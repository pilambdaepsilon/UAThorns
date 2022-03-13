/********************************
 * CONVERT ET ID TO IllinoisGRMHD
 *
 * Written in 2014 by Zachariah B. Etienne
 *
 * Sets metric & MHD variables needed
 * by IllinoisGRMHD, converting from
 * HydroBase and ADMBase.
 ********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "IllinoisGRMHD_headers.h"

#define cactusM (5.028916268544129e-34)    /*  1/g  */
#define cactusL (6.772400341316594e-06)    /*  1/cm */
#define cactusT (2.0303145448833407e5)     /*  1/s  */
#define cactusV (1.0/(cactusL*cactusL*cactusL))
#define RHO0    (1.0e15)

extern "C" void set_IllinoisGRMHD_metric_GRMHD_variables_based_on_HydroBase_and_ADMBase_variables(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  printf("I am starting initial data\n");

  gamma_tab[0] = Gamma_Initial;
  k_tab[0]     = K_Initial;

  if(rho_b_atm > 1e199) {
    CCTK_VError(VERR_DEF_PARAMS, "You MUST set rho_b_atm to some reasonable value in your param.ccl file.\n");
  }

  // Convert ADM variables (from ADMBase) to the BSSN-based variables expected by this routine.
  IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij(cctkGH,cctk_lsh,  gxx,gxy,gxz,gyy,gyz,gzz,alp,
                                                                              gtxx,gtxy,gtxz,gtyy,gtyz,gtzz,
                                                                              gtupxx,gtupxy,gtupxz,gtupyy,gtupyz,gtupzz,
                                                                              phi_bssn,psi_bssn,lapm1);

  if(pure_hydro_run) {
#pragma omp parallel for
    for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
          int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
          Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)]=0;
          Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)]=0;
          Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)]=0;
          Aphi[index]=0;
        }
  }
  double BARYMASS = 1.659e-24;                            //Same as in Cook code params.f
  double GNEWT = 6.672e-8;                                //Same as in Cook code params.f
  double CSPEED = 2.99792458e10;                          //Same as in Cook code params.f
  double barydconv = BARYMASS * GNEWT/(CSPEED*CSPEED);    //Same as in Cook code RotNS.f (makeEOS routine)
  double enrgdconv = GNEWT/(CSPEED*CSPEED);               //Same as in Cook code RotNS.f (makeEOS routine)
  double presconv = GNEWT/(CSPEED*CSPEED*CSPEED*CSPEED);  //Same as in Cook code RotNS.f (makeEOS routine)
  double SCALEK = 1.0/(RHO0*enrgdconv); 

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
	int index=CCTK_GFINDEX3D(cctkGH,i,j,k);

	rho_b[index] = rho[index]; rho_b[index] = MAX(rho_b[index], 1e-10); 
	P[index] = press[index];
        if(init_real){
           Temp[index] = 0.1;     //just initializing, we get temperature from other quantities later
           Wlorentz[index] = 1.0; //just initializes... this is fixed after C2P
           Ent[index] = 1.0;
           prim_eps[index]  = eps[index];
           Ye[index]   = Y_e[index];Ye[index] = MAX(Ye[index], 0.036);Ye[index] = MIN(Ye[index], 0.25);
           const int keytemp = 0; //assume temperature is unknown
	   const int eoskey = 4;
           const double precEOS = *c2p_eos_prec;
           double xeps = prim_eps[index];
           double xye = Ye[index];
           double xprs = P[index];
           double xrho = rho_b[index];
           double xtemp = Temp[index];
	   double P_cold = 0.0;
           int anyerr = 0;
           int keyerr = 0;
           const int npoints = 1;
	   double cH = ((xeps + xprs/xrho) + 1.0)*xrho; //given: rho, eps, P (and as such enthalpy), and Ye... get: T
           EOS_Omni_press_from_rhoenthalpy(eoskey, keytemp, precEOS, 1, &xrho, &xeps, &xtemp, &xye, &xprs, &cH, &anyerr, &keyerr); //get new T using this cH (also gives consistent P)
	   P[index] = xprs;                             //reset P to be consistent with cH(eps), rho, Ye
           Temp[index] = xtemp;                         
           if (rho_b[index] <=  *c2p_eos_rho_atmo){
             P[index]          = *c2p_eos_P_atmo;
             Temp[index]  = *c2p_eos_T_atmo;
             Ye[index]    = *c2p_eos_ye_atmo;
             prim_eps[index]   = *c2p_eos_eps_atmo;
             rho_b[index]      = *c2p_eos_rho_atmo;
             Ent[index]      = *c2p_eos_rho_atmo;
           }
        }
        else{
           Temp[index] = *c2p_eos_T_atmo;
           prim_eps[index]  = P[index]/rho_b[index];
           Ye[index]   = rho_b[index];
           if (rho_b[index] <=  *c2p_eos_rho_atmo){
             P[index]          = *c2p_eos_P_atmo;
             Temp[index]  = *c2p_eos_T_atmo;
             Ye[index]    = *c2p_eos_rho_atmo;
             prim_eps[index]   = *c2p_eos_eps_atmo;
             rho_b[index]      = *c2p_eos_rho_atmo;
           }

        }

        if (*c2p_eos_eoskey == 2){
          // P = (\Gamma - 1) rho epsilon
          // -> \Gamma = P/(rho epsilon) + 1
          double measured_gamma = ( P[index]/(rho_b[index] * prim_eps[index]) + 1.0 );
          if(rho_b[index]>rho_b_atm && fabs(gamma_th - measured_gamma)/gamma_th > 1e-2)
            CCTK_VError(VERR_DEF_PARAMS,
                        "Expected simple gamma law with gamma_th=%.15e, but found a point with gamma law such that gamma_th=%.15e. error = %e| rb=%e rbatm=%e P=%e, index=%d\n",
                        gamma_th,measured_gamma,(gamma_th-measured_gamma)/gamma_th,rho_b[index],rho_b_atm,P[index], index);
        }

	Ax[index] = Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)];
	Ay[index] = Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)];
	Az[index] = Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)];
	psi6phi[index] = Aphi[index];

	double ETvx = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)];
	double ETvy = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)];
	double ETvz = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)];

	// IllinoisGRMHD defines v^i = u^i/u^0.

	// Meanwhile, the ET/HydroBase formalism, called the Valencia
	// formalism, splits the 4 velocity into a purely spatial part
	// and a part that is normal to the spatial hypersurface:
	// u^a = G (n^a + U^a), (Eq. 14 of arXiv:1304.5544; G=W, U^a=v^a)
	// where n^a is the unit normal vector to the spatial hypersurface,
	// n_a = {-\alpha,0,0,0}, and U^a is the purely spatial part, which
	// is defined in HydroBase as the vel[] vector gridfunction.
	// Then u^a n_a = - \alpha u^0 = G n^a n_a = -G, and
	// of course \alpha u^0 = 1/sqrt(1+γ^ij u_i u_j) = \Gamma,
	// the standard Lorentz factor.

	// Note that n^i = - \beta^i / \alpha, so
	// u^a = \Gamma (n^a + U^a)
	// -> u^i = \Gamma ( U^i - \beta^i / \alpha )
	// which implies
	// v^i = u^i/u^0
	//     = \Gamma/u^0 ( U^i - \beta^i / \alpha ) <- \Gamma = \alpha u^0
	//     = \alpha ( U^i - \beta^i / \alpha )
	//     = \alpha U^i - \beta^i

	vx[index] = alp[index]*ETvx - betax[index];
	vy[index] = alp[index]*ETvy - betay[index];
	vz[index] = alp[index]*ETvz - betaz[index];

      }

  // Neat feature for debugging: Add a roundoff-error perturbation
  //    to the initial data.
  // Set random_pert variable to ~1e-14 for a random 15th digit
  //    perturbation.
  if (random_pert != 0){
    srandom(random_seed);
    if (pert_only_rho){
      for(int k=0;k<cctk_lsh[2];k++)
        for(int j=0;j<cctk_lsh[1];j++)
          for(int i=0;i<cctk_lsh[0];i++) {
            int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
            double pert = (random_pert*(double)rand() / RAND_MAX);
            double one_plus_pert=(1.0+pert);
            rho[index]*=one_plus_pert;
          }
    }else{
      for(int k=0;k<cctk_lsh[2];k++)
        for(int j=0;j<cctk_lsh[1];j++)
          for(int i=0;i<cctk_lsh[0];i++) {
            int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
            double pert = (random_pert*(double)rand() / RAND_MAX);
            double one_plus_pert=(1.0+pert);
            rho[index]*=one_plus_pert;
            vx[index]*=one_plus_pert;
            vy[index]*=one_plus_pert;
            vz[index]*=one_plus_pert;
            Ye[index]*=one_plus_pert;
            prim_eps[index]*=one_plus_pert;
            Temp[index]*=one_plus_pert;

            psi6phi[index]*=one_plus_pert;
            Ax[index]*=one_plus_pert;
            Ay[index]*=one_plus_pert;
            Az[index]*=one_plus_pert;
          }
    }
  }

  // Next compute B & B_stagger from A_i. Note that this routine also depends on
  //   the psi_bssn[] gridfunction being set to exp(phi).

  double dxi = 1.0/CCTK_DELTA_SPACE(0);
  double dyi = 1.0/CCTK_DELTA_SPACE(1);
  double dzi = 1.0/CCTK_DELTA_SPACE(2);

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {
        // Look Mom, no if() statements!
        int shiftedim1 = (i-1)*(i!=0); // This way, i=0 yields shiftedim1=0 and shiftedi=1, used below for our COPY boundary condition.
        int shiftedi   = shiftedim1+1;

        int shiftedjm1 = (j-1)*(j!=0);
        int shiftedj   = shiftedjm1+1;

        int shiftedkm1 = (k-1)*(k!=0);
        int shiftedk   = shiftedkm1+1;

        int index,indexim1,indexjm1,indexkm1;

        int actual_index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        double Psi = psi_bssn[actual_index];
        double Psim3 = 1.0/(Psi*Psi*Psi);

        // For the lower boundaries, the following applies a "copy"
        //    boundary condition on Bi_stagger where needed.
        //    E.g., Bx_stagger(i,jmin,k) = Bx_stagger(i,jmin+1,k)
        //    We find the copy BC works better than extrapolation.
        // For the upper boundaries, we do the following copy:
        //    E.g., Psi(imax+1,j,k)=Psi(imax,j,k)
        /**************/
        /* Bx_stagger */
        /**************/

        index    = CCTK_GFINDEX3D(cctkGH,i,shiftedj,shiftedk);
        indexjm1 = CCTK_GFINDEX3D(cctkGH,i,shiftedjm1,shiftedk);
        indexkm1 = CCTK_GFINDEX3D(cctkGH,i,shiftedj,shiftedkm1);
        // Set Bx_stagger = \partial_y A_z - partial_z A_y
        // "Grid" Ax(i,j,k) is actually Ax(i,j+1/2,k+1/2)
        // "Grid" Ay(i,j,k) is actually Ay(i+1/2,j,k+1/2)
        // "Grid" Az(i,j,k) is actually Ay(i+1/2,j+1/2,k)
        // Therefore, the 2nd order derivative \partial_z A_y at (i+1/2,j,k) is:
        //          ["Grid" Ay(i,j,k) - "Grid" Ay(i,j,k-1)]/dZ
        Bx_stagger[actual_index] = (Az[index]-Az[indexjm1])*dyi - (Ay[index]-Ay[indexkm1])*dzi;

        // Now multiply Bx and Bx_stagger by 1/sqrt(gamma(i+1/2,j,k)]) = 1/sqrt(1/2 [gamma + gamma_ip1]) = exp(-6 x 1/2 [phi + phi_ip1] )
        int imax_minus_i = (cctk_lsh[0]-1)-i;
        int indexip1jk = CCTK_GFINDEX3D(cctkGH,i + ( (imax_minus_i > 0) - (0 > imax_minus_i) ),j,k);
        double Psi_ip1 = psi_bssn[indexip1jk];
        Bx_stagger[actual_index] *= Psim3/(Psi_ip1*Psi_ip1*Psi_ip1);

        /**************/
        /* By_stagger */
        /**************/

        index    = CCTK_GFINDEX3D(cctkGH,shiftedi,j,shiftedk);
        indexim1 = CCTK_GFINDEX3D(cctkGH,shiftedim1,j,shiftedk);
        indexkm1 = CCTK_GFINDEX3D(cctkGH,shiftedi,j,shiftedkm1);
        // Set By_stagger = \partial_z A_x - \partial_x A_z
        By_stagger[actual_index] = (Ax[index]-Ax[indexkm1])*dzi - (Az[index]-Az[indexim1])*dxi;

        // Now multiply By and By_stagger by 1/sqrt(gamma(i,j+1/2,k)]) = 1/sqrt(1/2 [gamma + gamma_jp1]) = exp(-6 x 1/2 [phi + phi_jp1] )
        int jmax_minus_j = (cctk_lsh[1]-1)-j;
        int indexijp1k = CCTK_GFINDEX3D(cctkGH,i,j + ( (jmax_minus_j > 0) - (0 > jmax_minus_j) ),k);
        double Psi_jp1 = psi_bssn[indexijp1k];
        By_stagger[actual_index] *= Psim3/(Psi_jp1*Psi_jp1*Psi_jp1);


        /**************/
        /* Bz_stagger */
        /**************/

        index    = CCTK_GFINDEX3D(cctkGH,shiftedi,shiftedj,k);
        indexim1 = CCTK_GFINDEX3D(cctkGH,shiftedim1,shiftedj,k);
        indexjm1 = CCTK_GFINDEX3D(cctkGH,shiftedi,shiftedjm1,k);
        // Set Bz_stagger = \partial_x A_y - \partial_y A_x
        Bz_stagger[actual_index] = (Ay[index]-Ay[indexim1])*dxi - (Ax[index]-Ax[indexjm1])*dyi;

        // Now multiply Bz_stagger by 1/sqrt(gamma(i,j,k+1/2)]) = 1/sqrt(1/2 [gamma + gamma_kp1]) = exp(-6 x 1/2 [phi + phi_kp1] )
        int kmax_minus_k = (cctk_lsh[2]-1)-k;
        int indexijkp1 = CCTK_GFINDEX3D(cctkGH,i,j,k + ( (kmax_minus_k > 0) - (0 > kmax_minus_k) ));
        double Psi_kp1 = psi_bssn[indexijkp1];
        Bz_stagger[actual_index] *= Psim3/(Psi_kp1*Psi_kp1*Psi_kp1);

      }

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {
        // Look Mom, no if() statements!
        int shiftedim1 = (i-1)*(i!=0); // This way, i=0 yields shiftedim1=0 and shiftedi=1, used below for our COPY boundary condition.
        int shiftedi   = shiftedim1+1;

        int shiftedjm1 = (j-1)*(j!=0);
        int shiftedj   = shiftedjm1+1;

        int shiftedkm1 = (k-1)*(k!=0);
        int shiftedk   = shiftedkm1+1;

        int index,indexim1,indexjm1,indexkm1;

        int actual_index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // For the lower boundaries, the following applies a "copy"
        //    boundary condition on Bi and Bi_stagger where needed.
        //    E.g., Bx(imin,j,k) = Bx(imin+1,j,k)
        //    We find the copy BC works better than extrapolation.
        /******/
        /* Bx */
        /******/
        index = CCTK_GFINDEX3D(cctkGH,shiftedi,j,k);
        indexim1 = CCTK_GFINDEX3D(cctkGH,shiftedim1,j,k);
        // Set Bx = 0.5 ( Bx_stagger + Bx_stagger_im1 )
        // "Grid" Bx_stagger(i,j,k) is actually Bx_stagger(i+1/2,j,k)
        Bx[actual_index] = 0.5 * ( Bx_stagger[index] + Bx_stagger[indexim1] );

        /******/
        /* By */
        /******/
        index = CCTK_GFINDEX3D(cctkGH,i,shiftedj,k);
        indexjm1 = CCTK_GFINDEX3D(cctkGH,i,shiftedjm1,k);
        // Set By = 0.5 ( By_stagger + By_stagger_im1 )
        // "Grid" By_stagger(i,j,k) is actually By_stagger(i,j+1/2,k)
        By[actual_index] = 0.5 * ( By_stagger[index] + By_stagger[indexjm1] );

        /******/
        /* Bz */
        /******/
        index = CCTK_GFINDEX3D(cctkGH,i,j,shiftedk);
        indexkm1 = CCTK_GFINDEX3D(cctkGH,i,j,shiftedkm1);
        // Set Bz = 0.5 ( Bz_stagger + Bz_stagger_im1 )
        // "Grid" Bz_stagger(i,j,k) is actually Bz_stagger(i,j+1/2,k)
        Bz[actual_index] = 0.5 * ( Bz_stagger[index] + Bz_stagger[indexkm1] );
      }

  // FIXME: IllinoisGRMHD's Conservative-to-Primitive solver
  //   (a.k.a., C2P or con2prim) only implements single gamma-law EOS.
  // Also, not compatible with EOS driver in ET, so EOS parameters must
  //   be specified for both initial data thorns AND IllinoisGRMHD
  // TODO: Incorporate checks to ensure compatibility with ID.
  //   Alternatively, read in EOS stuff from an ET EOS driver.
  eos_struct eos;
  eos.neos=neos;
  eos.K_poly=K_poly;
  eos.rho_tab[0]=rho_tab[0];
  eos.P_tab[0]=P_tab[0];
  eos.gamma_th=gamma_th;
  eos.eps_tab[0]=eps_tab[0];
  eos.k_tab[0]=k_tab[0];   eos.k_tab[1]=k_tab[1];
  eos.gamma_tab[0]=gamma_tab[0]; eos.gamma_tab[1]=gamma_tab[1];

  // Finally, enforce limits on primitives & compute conservative variables.
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {
	static const int zero_int=0;
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	int ww;

        double PRIMS[MAXNUMVARS];
        ww=0;
        PRIMS[ww]  = rho_b[index]; ww++;
        PRIMS[ww]  = P[index];     ww++;
        PRIMS[ww]  = vx[index];    ww++;
        PRIMS[ww]  = vy[index];    ww++;
        PRIMS[ww]  = vz[index];    ww++;
        PRIMS[ww]  = Bx[index];    ww++;
        PRIMS[ww]  = By[index];    ww++;
        PRIMS[ww]  = Bz[index];
        ww = 17;
        PRIMS[ww]  = Ye[index];    ww++;
        PRIMS[ww]  = prim_eps[index];   ww++;
        PRIMS[ww]  = Temp[index];  ww++;
        PRIMS[ww]  = Wlorentz[index];  ww++;
        PRIMS[ww]  = Ent[index];  ww++;

	double METRIC[NUMVARS_FOR_METRIC],dummy=0;
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

	double CONSERVS[NUM_CONSERVS] = {0,0,0,0,0,0};
	double g4dn[4][4];
	double g4up[4][4];
	double TUPMUNU[10],TDNMUNU[10];

        struct output_stats stats; stats.failure_checker=0;
	IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(zero_int,PRIMS,stats,eos,METRIC,g4dn,g4up,TUPMUNU,TDNMUNU,CONSERVS,index);

	rho_b[index] = PRIMS[RHOB];
	Ye[index]    = PRIMS[YE];
        prim_eps[index] = PRIMS[EPS];
        Temp[index] = PRIMS[TEMP];
	P[index]     = PRIMS[PRESSURE];
	vx[index]    = PRIMS[VX];
	vy[index]    = PRIMS[VY];
	vz[index]    = PRIMS[VZ];

	rho_star[index] = CONSERVS[RHOSTAR];
	Yet[index]      = CONSERVS[YET];
	mhd_st_x[index] = CONSERVS[STILDEX];
	mhd_st_y[index] = CONSERVS[STILDEY];
	mhd_st_z[index] = CONSERVS[STILDEZ];
	tau[index]      = CONSERVS[TAUENERGY];

        if(update_Tmunu) {
          ww=0;
          eTtt[index] += TDNMUNU[ww]; ww++;
          eTtx[index] += TDNMUNU[ww]; ww++;
          eTty[index] += TDNMUNU[ww]; ww++;
          eTtz[index] += TDNMUNU[ww]; ww++;
          eTxx[index] += TDNMUNU[ww]; ww++;
          eTxy[index] += TDNMUNU[ww]; ww++;
          eTxz[index] += TDNMUNU[ww]; ww++;
          eTyy[index] += TDNMUNU[ww]; ww++;
          eTyz[index] += TDNMUNU[ww]; ww++;
          eTzz[index] += TDNMUNU[ww];
        }
      }

  printf("I am done with initial data\n");

  }
