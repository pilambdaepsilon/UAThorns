/* We evolve forward in time a set of functions called the
 *  "conservative variables", and any time the conserv's
 *  are updated, we must solve for the primitive variables
 *  (rho, pressure, velocities) using a Newton-Raphson
 *  technique, before reconstructing & evaluating the RHSs
 *  of the MHD equations again.
 *
 * This file contains the driver routine for this Newton-
 *  Raphson solver. Truncation errors in conservative
 *  variables can lead to no physical solutions in
 *  primitive variables. We correct for these errors here
 *  through a number of tricks described in the appendices
 *  of http://arxiv.org/pdf/1112.0568.pdf.
 *
 * This is a wrapper for the 2d solver of Noble et al. See
 *  harm_utoprim_2d.c for references and copyright notice
 *  for that solver. This wrapper was primarily written by
 *  Zachariah Etienne & Yuk Tung Liu, in 2011-2013.
 *
 * For optimal compatibility, this wrapper is licensed under
 *  the GPL v2 or any later version.
 *
 * Note that this code assumes a simple gamma law for the
 *  moment, though it would be easy to extend to a piecewise
 *  polytrope. */


#include "cctk.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/time.h>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#include "IllinoisGRMHD_headers.h"
#include "harm_primitives_headers.h"
#include "inlined_functions.C"
#include "apply_tau_floor__enforce_limits_on_primitives_and_recompute_conservs.C"

extern "C" void IllinoisGRMHD_conserv_to_prims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // We use proper C++ here, for file I/O later.
  using namespace std;

  // FIXME: only for single gamma-law EOS.
  eos_struct eos;
  eos.neos=neos;
  eos.K_poly=K_poly;
  eos.rho_tab[0]=rho_tab[0];
  eos.P_tab[0]=P_tab[0];
  eos.gamma_th=gamma_th;
  eos.eps_tab[0]=eps_tab[0];
  eos.k_tab[0]=k_tab[0];   eos.k_tab[1]=k_tab[1];
  eos.gamma_tab[0]=gamma_tab[0]; eos.gamma_tab[1]=gamma_tab[1];

  // These BSSN-based variables are not evolved, and so are not defined anywhere that the grid has moved.
  // Here we convert ADM variables (from ADMBase) to the BSSN-based variables expected by this routine.
  IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij(cctkGH,cctk_lsh,  gxx,gxy,gxz,gyy,gyz,gzz,alp,
                                                                gtxx,gtxy,gtxz,gtyy,gtyz,gtzz,
                                                                gtupxx,gtupxy,gtupxz,gtupyy,gtupyz,gtupzz,
                                                                phi_bssn,psi_bssn,lapm1);

  if(CCTK_EQUALS(Symmetry,"equatorial")) {
    // SET SYMMETRY GHOSTZONES ON ALL CONSERVATIVE VARIABLES!
    int ierr=0;
    ierr+=CartSymGN(cctkGH,"IllinoisGRMHD::grmhd_conservatives");
    // FIXME: UGLY. Filling metric ghostzones is needed for, e.g., Cowling runs.
    ierr+=CartSymGN(cctkGH,"lapse::lapse_vars");
    ierr+=CartSymGN(cctkGH,"bssn::BSSN_vars");
    ierr+=CartSymGN(cctkGH,"bssn::BSSN_AH");
    ierr+=CartSymGN(cctkGH,"shift::shift_vars");
    if(ierr!=0) CCTK_VError(VERR_DEF_PARAMS,"IllinoisGRMHD ERROR (grep for it, foo!)  :(");
  }

  //Start the timer, so we can benchmark the primitives solver during evolution.
  //  Slower solver -> harder to find roots -> things may be going crazy!
  //FIXME: Replace this timing benchmark with something more meaningful, like the avg # of Newton-Raphson iterations per gridpoint!
  /*
    struct timeval start, end;
    long mtime, seconds, useconds;
    gettimeofday(&start, NULL);
  */

  int gamma_equals2 = 1;
  if (fabs(gamma_th-2.0) > 1.e-10) gamma_equals2 = 0;

  int failures=0,font_fixes=0,vel_limited_ptcount=0;
  int pointcount=0;
  int failures_inhoriz=0;
  int pointcount_inhoriz=0;

  int pressure_cap_hit=0;

  CCTK_REAL error_int_numer=0,error_int_denom=0;

  int imin=0,jmin=0,kmin=0;
  int imax=cctk_lsh[0],jmax=cctk_lsh[1],kmax=cctk_lsh[2];

  int rho_star_fix_applied=0;
  long n_iter=0;
#pragma omp parallel for reduction(+:failures,vel_limited_ptcount,font_fixes,pointcount,failures_inhoriz,pointcount_inhoriz,error_int_numer,error_int_denom,pressure_cap_hit,rho_star_fix_applied,n_iter) schedule(static)
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        int ww;
        CCTK_REAL METRIC[NUMVARS_FOR_METRIC],dummy=0;
        ww=0;
        // FIXME: NECESSARY?
        //psi_bssn[index] = exp(phi[index]);
        METRIC[ww] = phi_bssn[index];ww++;
        METRIC[ww] = dummy;          ww++; // Don't need to set psi.
        METRIC[ww] = gtxx[index];    ww++;
        METRIC[ww] = gtxy[index];    ww++;
        METRIC[ww] = gtxz[index];    ww++;
        METRIC[ww] = gtyy[index];    ww++;
        METRIC[ww] = gtyz[index];    ww++;
        METRIC[ww] = gtzz[index];    ww++;
        METRIC[ww] = lapm1[index];   ww++;
        METRIC[ww] = betax[index];   ww++;
        METRIC[ww] = betay[index];   ww++;
        METRIC[ww] = betaz[index];   ww++;
        METRIC[ww] = gtupxx[index];  ww++;
        METRIC[ww] = gtupyy[index];  ww++;
        METRIC[ww] = gtupzz[index];  ww++;
        METRIC[ww] = gtupxy[index];  ww++;
        METRIC[ww] = gtupxz[index];  ww++;
        METRIC[ww] = gtupyz[index];  ww++;

        CCTK_REAL PRIMS[MAXNUMVARS];
        ww=0;
        PRIMS[ww] = rho_b[index]; ww++;
        PRIMS[ww] = P[index];     ww++;
        PRIMS[ww] = vx[index];    ww++;
        PRIMS[ww] = vy[index];    ww++;
        PRIMS[ww] = vz[index];    ww++;
        PRIMS[ww] = Bx[index];    ww++;
        PRIMS[ww] = By[index];    ww++;
        PRIMS[ww] = Bz[index];    ww++;
        PRIMS[YE] = Ye[index];    ww++;
        PRIMS[EPS] = prim_eps[index]; ww++;
        PRIMS[TEMP] = Temp[index]; ww++;
        PRIMS[WLORENTZ] = Wlorentz[index]; ww++;
        PRIMS[ENTROPY] = Ent[index]; ww++;

        CCTK_REAL CONSERVS[NUM_CONSERVS] = {rho_star[index], mhd_st_x[index],mhd_st_y[index],mhd_st_z[index],tau[index],Yet[index]};

        CCTK_REAL METRIC_LAP_PSI4[NUMVARS_METRIC_AUX];
        SET_LAPSE_PSI4(METRIC_LAP_PSI4,METRIC);

        CCTK_REAL METRIC_PHYS[NUMVARS_FOR_METRIC];
        METRIC_PHYS[GXX]   = METRIC[GXX]*METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GXY]   = METRIC[GXY]*METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GXZ]   = METRIC[GXZ]*METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GYY]   = METRIC[GYY]*METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GYZ]   = METRIC[GYZ]*METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GZZ]   = METRIC[GZZ]*METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GUPXX] = METRIC[GUPXX]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPXY] = METRIC[GUPXY]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPXZ] = METRIC[GUPXZ]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPYY] = METRIC[GUPYY]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPYZ] = METRIC[GUPYZ]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPZZ] = METRIC[GUPZZ]*METRIC_LAP_PSI4[PSIM4];

        CCTK_REAL TUPMUNU[10],TDNMUNU[10];

        CCTK_REAL shift_xL = METRIC_PHYS[GXX]*METRIC[SHIFTX] + METRIC_PHYS[GXY]*METRIC[SHIFTY] + METRIC_PHYS[GXZ]*METRIC[SHIFTZ];
        CCTK_REAL shift_yL = METRIC_PHYS[GXY]*METRIC[SHIFTX] + METRIC_PHYS[GYY]*METRIC[SHIFTY] + METRIC_PHYS[GYZ]*METRIC[SHIFTZ];
        CCTK_REAL shift_zL = METRIC_PHYS[GXZ]*METRIC[SHIFTX] + METRIC_PHYS[GYZ]*METRIC[SHIFTY] + METRIC_PHYS[GZZ]*METRIC[SHIFTZ];
        CCTK_REAL beta2L   = shift_xL*METRIC[SHIFTX] + shift_yL*METRIC[SHIFTY] + shift_zL*METRIC[SHIFTZ];

        // Compute 4-metric, both g_{\mu \nu} and g^{\mu \nu}.
        // This is for computing T_{\mu \nu} and T^{\mu \nu}. Also the HARM con2prim lowlevel function requires them.
        CCTK_REAL g4dn[4][4],g4up[4][4];
        g4dn[0][0] = -SQR(METRIC_LAP_PSI4[LAPSE]) + beta2L;
        g4dn[0][1] = g4dn[1][0] = shift_xL;
        g4dn[0][2] = g4dn[2][0] = shift_yL;
        g4dn[0][3] = g4dn[3][0] = shift_zL;
        g4dn[1][1]              = METRIC_PHYS[GXX];
        g4dn[1][2] = g4dn[2][1] = METRIC_PHYS[GXY];
        g4dn[1][3] = g4dn[3][1] = METRIC_PHYS[GXZ];
        g4dn[2][2]              = METRIC_PHYS[GYY];
        g4dn[2][3] = g4dn[3][2] = METRIC_PHYS[GYZ];
        g4dn[3][3]              = METRIC_PHYS[GZZ];

        CCTK_REAL alpha_inv_squared=SQR(METRIC_LAP_PSI4[LAPSEINV]);
        g4up[0][0] = -1.0*alpha_inv_squared;
        g4up[0][1] = g4up[1][0] = METRIC[SHIFTX]*alpha_inv_squared;
        g4up[0][2] = g4up[2][0] = METRIC[SHIFTY]*alpha_inv_squared;
        g4up[0][3] = g4up[3][0] = METRIC[SHIFTZ]*alpha_inv_squared;
        g4up[1][1]              = METRIC_PHYS[GUPXX] - METRIC[SHIFTX]*METRIC[SHIFTX]*alpha_inv_squared;
        g4up[1][2] = g4up[2][1] = METRIC_PHYS[GUPXY] - METRIC[SHIFTX]*METRIC[SHIFTY]*alpha_inv_squared;
        g4up[1][3] = g4up[3][1] = METRIC_PHYS[GUPXZ] - METRIC[SHIFTX]*METRIC[SHIFTZ]*alpha_inv_squared;
        g4up[2][2]              = METRIC_PHYS[GUPYY] - METRIC[SHIFTY]*METRIC[SHIFTY]*alpha_inv_squared;
        g4up[2][3] = g4up[3][2] = METRIC_PHYS[GUPYZ] - METRIC[SHIFTY]*METRIC[SHIFTZ]*alpha_inv_squared;
        g4up[3][3]              = METRIC_PHYS[GUPZZ] - METRIC[SHIFTZ]*METRIC[SHIFTZ]*alpha_inv_squared;

        //FIXME: might slow down the code.
        if(isnan(CONSERVS[RHOSTAR]*CONSERVS[STILDEX]*CONSERVS[STILDEY]*CONSERVS[STILDEZ]*CONSERVS[TAUENERGY]*CONSERVS[YET]*PRIMS[BX_CENTER]*PRIMS[BY_CENTER]*PRIMS[BZ_CENTER])) {
          CCTK_VInfo(CCTK_THORNSTRING,"NAN FOUND: i,j,k = %d %d %d, x,y,z = %e %e %e , index=%d st_i = %e %e %e, rhostar = %e, Yet = %e, tau = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",
                     i,j,k,x[index],y[index],z[index],index,
                     CONSERVS[STILDEX],CONSERVS[STILDEY],CONSERVS[STILDEZ],CONSERVS[RHOSTAR],CONSERVS[YET],CONSERVS[TAUENERGY],
                     PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER],METRIC_PHYS[GXX],METRIC_PHYS[GXY],METRIC_PHYS[GXZ],METRIC_PHYS[GYY],METRIC_PHYS[GYZ],METRIC_PHYS[GZZ],METRIC_LAP_PSI4[PSI6]);
        }

        // Here we use _flux variables as temp storage for original values of conservative variables.. This is used for debugging purposes only.
        rho_star_flux[index] = CONSERVS[RHOSTAR];
        st_x_flux[index]     = CONSERVS[STILDEX];
        st_y_flux[index]     = CONSERVS[STILDEY];
        st_z_flux[index]     = CONSERVS[STILDEZ];
        tau_flux[index]      = CONSERVS[TAUENERGY];
        Yet_flux[index]      = CONSERVS[YET];

        CCTK_REAL rho_star_orig = CONSERVS[RHOSTAR];
        CCTK_REAL mhd_st_x_orig = CONSERVS[STILDEX];
        CCTK_REAL mhd_st_y_orig = CONSERVS[STILDEY];
        CCTK_REAL mhd_st_z_orig = CONSERVS[STILDEZ];
        CCTK_REAL tau_orig      = CONSERVS[TAUENERGY];
        CCTK_REAL Yet_orig      = CONSERVS[YET];

        int check=0;
        struct output_stats stats;
        stats.n_iter=0;
        stats.vel_limited=0;
        stats.failure_checker=0;

        CCTK_REAL xco=x[index]; CCTK_REAL yco=y[index]; CCTK_REAL zco=z[index];
        if(CONSERVS[RHOSTAR]>0.0) {
          // Apply the tau floor
          apply_tau_floor(index,tau_atm,rho_b_atm,Psi6threshold,PRIMS,METRIC,METRIC_PHYS,METRIC_LAP_PSI4,stats,eos,  CONSERVS);

          stats.font_fixed=0;
          for(int ii=0;ii<3;ii++) {
            //USE OLD SOLVERS
            if(!use_ConservativeToPrimitive){
              check = harm_primitives_gammalaw_lowlevel(index,i,j,k,x,y,z,METRIC,METRIC_PHYS,METRIC_LAP_PSI4,
                                                      CONSERVS,PRIMS,  g4dn,g4up,   stats,eos);
              //better not to leave these new variables un-updated.
              PRIMS[YE] = CONSERVS[YET]/CONSERVS[RHOSTAR];
              PRIMS[EPS] = PRIMS[PRESSURE]/((eos.gamma_tab[0]-1.0)*PRIMS[RHOB]);
              PRIMS[TEMP] = 0.1;
            }else{      //USE NEW SOLVERS
              CCTK_REAL PRIMS_orig[MAXNUMVARS];
              ww=0;
              PRIMS_orig[ww] = rho_b[index]; ww++;
              PRIMS_orig[ww] = P[index];     ww++;
              PRIMS_orig[ww] = vx[index];    ww++;
              PRIMS_orig[ww] = vy[index];    ww++;
              PRIMS_orig[ww] = vz[index];    ww++;
              PRIMS_orig[ww] = Bx[index];    ww++;
              PRIMS_orig[ww] = By[index];    ww++;
              PRIMS_orig[ww] = Bz[index];    ww++;

              // IllinoisGRMHD defines conservatives with sqrtgamma inside, ConservativeToPrimitive does not
              CCTK_REAL sqrtgamma = METRIC_LAP_PSI4[PSI6];
              CCTK_REAL prim[14];
              CCTK_REAL con[10];
              con[C2P_D]      = CONSERVS[RHOSTAR] / sqrtgamma;
              con[C2P_S1_cov] = CONSERVS[STILDEX] / sqrtgamma;
              con[C2P_S2_cov] = CONSERVS[STILDEY] / sqrtgamma;
              con[C2P_S3_cov] = CONSERVS[STILDEZ] / sqrtgamma;
              con[C2P_TAU]    = CONSERVS[TAUENERGY] / sqrtgamma;
              con[C2P_B1_con] = PRIMS[BX_CENTER];
              con[C2P_B2_con] = PRIMS[BY_CENTER];
              con[C2P_B3_con] = PRIMS[BZ_CENTER];
              con[C2P_YE]     = CONSERVS[YET] / sqrtgamma;
              con[C2P_TEMP]   = Temp[index];
              // Initial guess
              prim[C2P_RHO]    = PRIMS[RHOB];
              prim[C2P_PRESS]  = PRIMS[PRESSURE];
              // Convert the velocities to Valencia
              CCTK_REAL one_over_alpha = 1./METRIC_LAP_PSI4[LAPSE];
              CCTK_REAL v1_con = one_over_alpha * (PRIMS[VX] + METRIC[SHIFTX]);
              CCTK_REAL v2_con = one_over_alpha * (PRIMS[VY] + METRIC[SHIFTY]);
              CCTK_REAL v3_con = one_over_alpha * (PRIMS[VZ] + METRIC[SHIFTZ]);
              // Lower index
              prim[C2P_v1_cov] = METRIC_PHYS[GXX] * v1_con + METRIC_PHYS[GXY] * v2_con + METRIC_PHYS[GXZ] * v3_con;
              prim[C2P_v2_cov] = METRIC_PHYS[GXY] * v1_con + METRIC_PHYS[GYY] * v2_con + METRIC_PHYS[GYZ] * v3_con;
              prim[C2P_v3_cov] = METRIC_PHYS[GXZ] * v1_con + METRIC_PHYS[GYZ] * v2_con + METRIC_PHYS[GZZ] * v3_con;

              prim[C2P_EPS]    = prim_eps[index];
              prim[C2P_B1_con] = PRIMS[BX_CENTER];
              prim[C2P_B2_con] = PRIMS[BY_CENTER];
              prim[C2P_B3_con] = PRIMS[BZ_CENTER];
              prim[C2P_YE]     = Ye[index];
              prim[C2P_TEMP]   = Temp[index];
              prim[C2P_WLORENTZ] = con[C2P_D] / PRIMS[RHOB];
              prim[C2P_ENT]    = 0.;
              prim[C2P_A_BAR]  = 0.;

              // Main routine
              check = Con2Prim(prim, *g4up, *g4dn, con, cctk_iteration, index,xco,yco,zco, METRIC, METRIC_PHYS, METRIC_LAP_PSI4);
              PRIMS[RHOB]      = prim[C2P_RHO];
              PRIMS[PRESSURE]  = prim[C2P_PRESS];
              PRIMS[EPS]       = prim[C2P_EPS];
              PRIMS[TEMP]      = prim[C2P_TEMP];
              PRIMS[WLORENTZ]  = prim[C2P_WLORENTZ];
              PRIMS[ENTROPY]   = prim[C2P_ENT];
              PRIMS[YE]        = CONSERVS[YET]/CONSERVS[RHOSTAR];
              //Raise index
              v1_con = METRIC_PHYS[GUPXX] * prim[C2P_v1_cov] + METRIC_PHYS[GUPXY] * prim[C2P_v2_cov] + METRIC_PHYS[GUPXZ] * prim[C2P_v3_cov];
              v2_con = METRIC_PHYS[GUPXY] * prim[C2P_v1_cov] + METRIC_PHYS[GUPYY] * prim[C2P_v2_cov] + METRIC_PHYS[GUPYZ] * prim[C2P_v3_cov];
              v3_con = METRIC_PHYS[GUPXZ] * prim[C2P_v1_cov] + METRIC_PHYS[GUPYZ] * prim[C2P_v2_cov] + METRIC_PHYS[GUPZZ] * prim[C2P_v3_cov];
              CCTK_REAL lap = METRIC_LAP_PSI4[LAPSE];
              // Convert to IllinoisGRMHD formulation
              PRIMS[VX]        = lap * v1_con - METRIC[SHIFTX];
              PRIMS[VY]        = lap * v2_con - METRIC[SHIFTY];
              PRIMS[VZ]        = lap * v3_con - METRIC[SHIFTZ];
            }

            if(check==0) ii=4;
            else stats.failure_checker+=100000;
          }
        } else {
          stats.failure_checker+=1;

          PRIMS[RHOB] = *c2p_eos_rho_atmo;
          PRIMS[TEMP] = *c2p_eos_T_atmo;
          PRIMS[YE]   = *c2p_eos_ye_atmo;
          PRIMS[EPS]  = *c2p_eos_eps_atmo;

          if(!use_ConservativeToPrimitive){
            // Set to atmosphere if rho_star<0.
            PRIMS[RHOB]      =rho_b_atm;
            if(gamma_equals2) {
              PRIMS[PRESSURE]=K_poly*SQR(rho_b_atm);
            } else {
              PRIMS[PRESSURE]=K_poly*pow(rho_b_atm,gamma_th);
            }
          }else{
            PRIMS[PRESSURE] = *c2p_eos_P_atmo;
//            PRIMS[RHOB]      =rho_b_atm;
            if(gamma_equals2) {
              PRIMS[PRESSURE]=K_poly*SQR(rho_b_atm);
            } else {
              PRIMS[PRESSURE]=K_poly*pow(rho_b_atm,gamma_th);
            }
          }
          PRIMS[VX]        =-METRIC[SHIFTX];
          PRIMS[VY]        =-METRIC[SHIFTY];
          PRIMS[VZ]        =-METRIC[SHIFTZ];

          rho_star_fix_applied++;
        }

        // Enforce limits on primitive variables and recompute conservatives.
        static const int already_computed_physical_metric_and_inverse=1;
        IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(already_computed_physical_metric_and_inverse,PRIMS,stats,eos,METRIC,g4dn,g4up, TUPMUNU,TDNMUNU,CONSERVS,index);

        rho_star[index] = CONSERVS[RHOSTAR];
        Yet[index]      = CONSERVS[YET];
        mhd_st_x[index] = CONSERVS[STILDEX];
        mhd_st_y[index] = CONSERVS[STILDEY];
        mhd_st_z[index] = CONSERVS[STILDEZ];
        tau[index] = CONSERVS[TAUENERGY];

        // Set primitives, and/or provide a better guess.
        rho_b[index]    = PRIMS[RHOB];
        Ye[index]       = PRIMS[YE];
        P[index]        = PRIMS[PRESSURE];
        vx[index]       = PRIMS[VX];
        vy[index]       = PRIMS[VY];
        vz[index]       = PRIMS[VZ];
        prim_eps[index] = PRIMS[EPS];
        Temp[index]     = PRIMS[TEMP];
        Wlorentz[index] = PRIMS[WLORENTZ];
        Ent[index]      = PRIMS[ENTROPY];

        if(update_Tmunu) {
          ww=0;
          eTtt[index] = TDNMUNU[ww]; ww++;
          eTtx[index] = TDNMUNU[ww]; ww++;
          eTty[index] = TDNMUNU[ww]; ww++;
          eTtz[index] = TDNMUNU[ww]; ww++;
          eTxx[index] = TDNMUNU[ww]; ww++;
          eTxy[index] = TDNMUNU[ww]; ww++;
          eTxz[index] = TDNMUNU[ww]; ww++;
          eTyy[index] = TDNMUNU[ww]; ww++;
          eTyz[index] = TDNMUNU[ww]; ww++;
          eTzz[index] = TDNMUNU[ww];
        }

        //Finally, we set h, the enthalpy:
        //CCTK_REAL eps = P[index]/rho_b[index]/(GAMMA-1.0);
        //h[index] = 1.0 + P[index]/rho_b[index] + eps;

        /***************************************************************************************************************************/
        // DIAGNOSTICS:
        //Pressure cap hit?
        /* FIXME
        CCTK_REAL P_cold = rho_b[index]*rho_b[index];
        if(P[index]/P_cold > 0.99*1e3 && rho_b[index]>100.0*rho_b_atm) {
          if(exp(phi[index]*6.0) <= Psi6threshold) pressure_cap_hit++;
        }
        */

        //Now we compute the difference between original & new conservatives, for diagnostic purposes:
        error_int_numer += fabs(tau[index] - tau_orig) + fabs(rho_star[index] - rho_star_orig) +
          fabs(mhd_st_x[index] - mhd_st_x_orig) + fabs(mhd_st_y[index] - mhd_st_y_orig) + fabs(mhd_st_z[index] - mhd_st_z_orig) + fabs(Yet[index] - Yet_orig);
        error_int_denom += tau_orig + rho_star_orig + fabs(mhd_st_x_orig) + fabs(mhd_st_y_orig) + fabs(mhd_st_z_orig) + fabs(Yet_orig);

        if(stats.font_fixed==1) font_fixes++;
        vel_limited_ptcount+=stats.vel_limited;
        if(check!=0) {
          failures++;
          if(exp(METRIC[PHI]*6.0)>Psi6threshold) {
            failures_inhoriz++;
            pointcount_inhoriz++;
          }
        }
        pointcount++;
        /***************************************************************************************************************************/
        failure_checker[index] = stats.failure_checker;
        n_iter += stats.n_iter;
      }

  /*
    gettimeofday(&end, NULL);

    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;

    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.999;  // We add 0.999 since mtime is a long int; this rounds up the result before setting the value.  Here, rounding down is incorrect.
    solutions per second: cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2] / ((CCTK_REAL)mtime/1000.0),
  */
  if(CCTK_Equals(verbose, "essential") || CCTK_Equals(verbose, "essential+iteration output")) {
    CCTK_VInfo(CCTK_THORNSTRING,"C2P: Lev: %d NumPts= %d | Fixes: Font= %d VL= %d rho*= %d | Failures: %d InHoriz= %d / %d | Error: %.3e, ErrDenom: %.3e | %.2f iters/gridpt | %1.3e patm",
               (int)GetRefinementLevel(cctkGH),
               pointcount,font_fixes,vel_limited_ptcount,rho_star_fix_applied,
               failures,
               failures_inhoriz,pointcount_inhoriz,
               error_int_numer/error_int_denom,error_int_denom,
               (double)n_iter/( (double)(cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]) ), c2p_eos_P_atmo);
  }

  if(pressure_cap_hit!=0) {
    //CCTK_VInfo(CCTK_THORNSTRING,"PRESSURE CAP HIT %d TIMES!  Outputting debug file!",pressure_cap_hit);
  }

  // Very useful con2prim debugger. If the primitives (con2prim) solver fails, this will output all data needed to
  //     debug where and why the solver failed. Strongly suggested for experimenting with new fixes.
  if(conserv_to_prims_debug==1) {

    ofstream myfile;
    char filename[100];
    srand(time(NULL));
    sprintf(filename,"primitives_debug-%e.dat",rho_star[CCTK_GFINDEX3D(cctkGH,3,15,6)]);
    //Alternative, for debugging purposes as well: sprintf(filename,"primitives_debug-%d.dat",rand());
    myfile.open (filename, ios::out | ios::binary);
    //myfile.open ("data.bin", ios::out | ios::binary);
    myfile.write((char*)cctk_lsh, 3*sizeof(int));

    myfile.write((char*)&rho_b_atm, 1*sizeof(CCTK_REAL));
    myfile.write((char*)&tau_atm, 1*sizeof(CCTK_REAL));

    myfile.write((char*)&Psi6threshold, 1*sizeof(CCTK_REAL));

    CCTK_REAL gamma_th=2.0;
    myfile.write((char*)&gamma_th, 1*sizeof(CCTK_REAL));
    myfile.write((char*)&neos,     1*sizeof(int));

    myfile.write((char*)gamma_tab, (neos+1)*sizeof(CCTK_REAL));
    myfile.write((char*)k_tab,     (neos+1)*sizeof(CCTK_REAL));

    myfile.write((char*)eps_tab,   neos*sizeof(CCTK_REAL));
    myfile.write((char*)rho_tab,   neos*sizeof(CCTK_REAL));
    myfile.write((char*)P_tab,     neos*sizeof(CCTK_REAL));

    int fullsize=cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
    myfile.write((char*)x,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)y,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)z,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)phi_bssn, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtxx, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtxy, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtxz, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtyy, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtyz, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtzz, (fullsize)*sizeof(CCTK_REAL));

    myfile.write((char*)gtupxx, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupxy, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupxz, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupyy, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupyz, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupzz, (fullsize)*sizeof(CCTK_REAL));

    myfile.write((char*)betax, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)betay, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)betaz, (fullsize)*sizeof(CCTK_REAL));

    myfile.write((char*)lapm1, (fullsize)*sizeof(CCTK_REAL));

    // HERE WE USE _flux variables as temp storage for original values of conservative variables.. This is used for debugging purposes only.
    myfile.write((char*)tau_flux,      (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)st_x_flux, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)st_y_flux, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)st_z_flux, (fullsize)*sizeof(CCTK_REAL));

    myfile.write((char*)rho_star_flux, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)Yet_flux, (fullsize)*sizeof(CCTK_REAL));

    myfile.write((char*)Bx,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)By,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)Bz,   (fullsize)*sizeof(CCTK_REAL));

    myfile.write((char*)vx,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)vy,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)vz,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)P,    (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)rho_b,(fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)Ye,   (fullsize)*sizeof(CCTK_REAL));

    int checker=1063; myfile.write((char*)&checker,sizeof(int));

    myfile.close();
    CCTK_VInfo(CCTK_THORNSTRING,"Finished writing...");
  }
}

#include "harm_primitives_lowlevel.C"
