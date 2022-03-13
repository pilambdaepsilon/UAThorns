/*****************************************
 * PPM Reconstruction Interface.
 * Zachariah B. Etienne (2013)
 *
 * This version of PPM implements the standard
 * Colella & Woodward PPM, though modified as in GRHydro
 * to have 3 ghostzones instead of 4.
 *****************************************/

#define MINUS2 0
#define MINUS1 1
#define PLUS0  2
#define PLUS1  3
#define PLUS2  4
#define MAXNUMINDICES 5
//    ^^^^^^^^^^^^^ Be _sure_ to define MAXNUMINDICES appropriately!

// You'll find the #define's for LOOP_DEFINE and SET_INDEX_ARRAYS inside:
#include "loop_defines_reconstruction.h"

static inline CCTK_REAL ftilde_compute(const int flux_dirn,CCTK_REAL U[MAXNUMVARS][MAXNUMINDICES]);
static inline CCTK_REAL slope_limit(CCTK_REAL dU,CCTK_REAL dUp1);
static inline void steepen_rho(const int eos_key,const int keytemp,const CCTK_REAL c2pprec,CCTK_REAL U[MAXNUMVARS][MAXNUMINDICES],CCTK_REAL slope_lim_dU[MAXNUMVARS][MAXNUMINDICES],
                               CCTK_REAL gamma_th,CCTK_REAL P_cold,CCTK_REAL gamma_cold,
                               CCTK_REAL *rho_br_ppm,CCTK_REAL *rho_bl_ppm);
static inline void compute_P_cold__gamma_cold(CCTK_REAL rho_b,eos_struct &eos,   CCTK_REAL &P_cold,CCTK_REAL &gamma_cold);
static inline void monotonize(CCTK_REAL U,CCTK_REAL &Ur,CCTK_REAL &Ul);

static void reconstruct_set_of_prims_PPM(const int eos_key,const int keytemp,const CCTK_REAL c2pprec,const cGH *cctkGH,const int *cctk_lsh,const int flux_dirn,const int num_prims_to_reconstruct,const int *which_prims_to_reconstruct,eos_struct &eos,
                                         gf_and_gz_struct *in_prims,gf_and_gz_struct *out_prims_r,gf_and_gz_struct *out_prims_l,CCTK_REAL *ftilde_gf, CCTK_REAL *temporary) {

  CCTK_REAL U[MAXNUMVARS][MAXNUMINDICES],dU[MAXNUMVARS][MAXNUMINDICES],slope_lim_dU[MAXNUMVARS][MAXNUMINDICES],
    Ur[MAXNUMVARS][MAXNUMINDICES],Ul[MAXNUMVARS][MAXNUMINDICES];
  int ijkgz_lo_hi[4][2];

  for(int ww=0;ww<num_prims_to_reconstruct;ww++) {
    int whichvar=which_prims_to_reconstruct[ww];

    if(in_prims[whichvar].gz_lo[flux_dirn]!=0 || in_prims[whichvar].gz_hi[flux_dirn]!=0) {
      CCTK_VError(VERR_DEF_PARAMS,"TOO MANY GZ'S! WHICHVAR=%d: %d %d %d : %d %d %d DIRECTION %d",whichvar,
		  in_prims[whichvar].gz_lo[1],in_prims[whichvar].gz_lo[2],in_prims[whichvar].gz_lo[3],
		  in_prims[whichvar].gz_hi[1],in_prims[whichvar].gz_hi[2],in_prims[whichvar].gz_hi[3],flux_dirn);
    }

    // *** LOOP 1: Interpolate to Ur and Ul, which are face values ***
    //  You will find that Ur depends on U at MINUS1,PLUS0, PLUS1,PLUS2, and
    //                     Ul depends on U at MINUS2,MINUS1,PLUS0,PLUS1.
    //  However, we define the below loop from MINUS2 to PLUS2. Why not split
    //     this up and get additional points? The reason is that later on,
    //     Ur and Ul depend on ftilde, which is defined from MINUS2 to PLUS2,
    //     so we would lose those points anyway.
    LOOP_DEFINE(2,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
      SET_INDEX_ARRAYS(-2,2,flux_dirn);
      /* *** LOOP 1a: READ INPUT *** */
      // Read in a primitive at all gridpoints between m = MINUS2 & PLUS2, where m's direction is given by flux_dirn. Store to U.
      for(int ii=MINUS2;ii<=PLUS2;ii++) U[whichvar][ii] = in_prims[whichvar].gf[index_arr[flux_dirn][ii]];

      /* *** LOOP 1b: DO COMPUTATION *** */
      /* First, compute simple dU = U(i) - U(i-1), where direction of i
       *         is given by flux_dirn, and U is a primitive variable:
       *         {rho_b,P,vx,vy,vz,Bx,By,Bz,Ye}. */
      // Note that for Ur and Ul at i, we must compute dU(i-1),dU(i),dU(i+1),
      //         and dU(i+2)
      dU[whichvar][MINUS1] = U[whichvar][MINUS1]- U[whichvar][MINUS2];
      dU[whichvar][PLUS0]  = U[whichvar][PLUS0] - U[whichvar][MINUS1];
      dU[whichvar][PLUS1]  = U[whichvar][PLUS1] - U[whichvar][PLUS0];
      dU[whichvar][PLUS2]  = U[whichvar][PLUS2] - U[whichvar][PLUS1];

      // Then, compute slope-limited dU, using MC slope limiter:
      slope_lim_dU[whichvar][MINUS1]=slope_limit(dU[whichvar][MINUS1],dU[whichvar][PLUS0]);
      slope_lim_dU[whichvar][PLUS0] =slope_limit(dU[whichvar][PLUS0], dU[whichvar][PLUS1]);
      slope_lim_dU[whichvar][PLUS1] =slope_limit(dU[whichvar][PLUS1], dU[whichvar][PLUS2]);

      // Finally, compute face values Ur and Ul based on the PPM prescription
      //   (Eq. A1 in http://arxiv.org/pdf/astro-ph/0503420.pdf, but using standard 1/6=(1.0/6.0) coefficient)
      // Ur[PLUS0] represents U(i+1/2)
      // We applied a simplification to the following line: Ur=U+0.5*(U(i+1)-U) + ... = 0.5*(U(i+1)+U) + ...
      Ur[whichvar][PLUS0] = 0.5*(U[whichvar][PLUS1] + U[whichvar][PLUS0] ) + (1.0/6.0)*(slope_lim_dU[whichvar][PLUS0]  - slope_lim_dU[whichvar][PLUS1]);
      // Ul[PLUS0] represents U(i-1/2)
      // We applied a simplification to the following line: Ul=U(i-1)+0.5*(U-U(i-1)) + ... = 0.5*(U+U(i-1)) + ...
      Ul[whichvar][PLUS0] = 0.5*(U[whichvar][PLUS0] + U[whichvar][MINUS1]) + (1.0/6.0)*(slope_lim_dU[whichvar][MINUS1] - slope_lim_dU[whichvar][PLUS0]);

      /* *** LOOP 1c: WRITE OUTPUT *** */
      // Store right face values to {rho_br,Pr,vxr,vyr,vzr,Bxr,Byr,Bzr},
      //    and left face values to {rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl}
      out_prims_r[whichvar].gf[index_arr[flux_dirn][PLUS0]] = Ur[whichvar][PLUS0];
      out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]] = Ul[whichvar][PLUS0];
    }

    // *** LOOP 2a: STEEPEN RHOB ***
    if(whichvar==RHOB) {
      LOOP_DEFINE(2,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
	SET_INDEX_ARRAYS(-2,2,flux_dirn);
	// Set rho and P separately, since within this loop,
	//  1) steepen_rho() depends on RHOB(MINUS2,MINUS1,PLUS0,PLUS1,PLUS2)

	// Read in all primitives between MINUS2 & PLUS2. Store to U.
	for(int ii=MINUS2;ii<=PLUS2;ii++) U[RHOB][ii]     = in_prims[RHOB    ].gf[index_arr[flux_dirn][ii]];
	for(int ii=MINUS2;ii<=PLUS2;ii++) U[YE][ii]       = in_prims[YE      ].gf[index_arr[flux_dirn][ii]];
	for(int ii=MINUS2;ii<=PLUS2;ii++) U[EPS][ii]      = in_prims[EPS     ].gf[index_arr[flux_dirn][ii]];
	for(int ii=MINUS2;ii<=PLUS2;ii++) U[TEMP][ii]     = in_prims[TEMP    ].gf[index_arr[flux_dirn][ii]];
	for(int ii=MINUS1;ii<=PLUS1;ii++) U[PRESSURE][ii] = in_prims[PRESSURE].gf[index_arr[flux_dirn][ii]];
	Ur[RHOB][PLUS0] = out_prims_r[RHOB].gf[index_arr[flux_dirn][PLUS0]];
	Ul[RHOB][PLUS0] = out_prims_l[RHOB].gf[index_arr[flux_dirn][PLUS0]];

	dU[whichvar][MINUS1] = U[whichvar][MINUS1]- U[whichvar][MINUS2];
	dU[whichvar][PLUS0]  = U[whichvar][PLUS0] - U[whichvar][MINUS1];
	dU[whichvar][PLUS1]  = U[whichvar][PLUS1] - U[whichvar][PLUS0];
	dU[whichvar][PLUS2]  = U[whichvar][PLUS2] - U[whichvar][PLUS1];

	slope_lim_dU[whichvar][MINUS1]=slope_limit(dU[whichvar][MINUS1],dU[whichvar][PLUS0]);
	slope_lim_dU[whichvar][PLUS1] =slope_limit(dU[whichvar][PLUS1], dU[whichvar][PLUS2]);

	// Steepen rho
	// DEPENDENCIES: RHOB face values, RHOB(MINUS2,MINUS1,PLUS0,PLUS1,PLUS2), P(MINUS1,PLUS0,PLUS1), and slope_lim_dU[RHOB](MINUS1,PLUS1)
	CCTK_REAL P_cold=0,gamma_cold=0;
	compute_P_cold__gamma_cold(U[RHOB][PLUS0],eos,  P_cold,gamma_cold);
	steepen_rho(eos_key,keytemp,c2pprec,U,slope_lim_dU,   eos.gamma_th,P_cold,gamma_cold,  Ur[RHOB],Ul[RHOB]);

	// Output rho
	out_prims_r[RHOB].gf[index_arr[flux_dirn][PLUS0]] = Ur[RHOB][PLUS0];
	out_prims_l[RHOB].gf[index_arr[flux_dirn][PLUS0]] = Ul[RHOB][PLUS0];
      }
    }
  }

  /* ORIGINAL PPM REQUIRES AT LEAST 4 GHOSTZONES, which can add
   *  significantly to the size of AMR ref. boundaries.
   *  To reduce to 3 ghostzones, we comment the following lines out:
   * if ((P[indexp1] - P[indexm1]) <= 0.0) {
   * f = MAX(ftilde,ftilde_p1);
   * } else {
   * f = MAX(ftilde,ftilde_m1);
   * }
   */

  // *** LOOP 3: FLATTEN BASED ON FTILDE AND MONOTONIZE ***
  for(int ww=0;ww<num_prims_to_reconstruct;ww++) {
    int whichvar=which_prims_to_reconstruct[ww];
    // ftilde() depends on P(MINUS2,MINUS1,PLUS1,PLUS2)
    LOOP_DEFINE(2,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
      SET_INDEX_ARRAYS(0,0,flux_dirn);

      U[whichvar][PLUS0]  = in_prims[whichvar].gf[index_arr[flux_dirn][PLUS0]];
      Ur[whichvar][PLUS0] = out_prims_r[whichvar].gf[index_arr[flux_dirn][PLUS0]];
      Ul[whichvar][PLUS0] = out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]];

      // ftilde_gf was computed in the function compute_ftilde_gf(), called before this routine
      CCTK_REAL ftilde = ftilde_gf[index_arr[flux_dirn][PLUS0]];
      // ...and then flatten (local operation)
      Ur[whichvar][PLUS0]   = U[whichvar][PLUS0]*ftilde + Ur[whichvar][PLUS0]*(1.0-ftilde);
      Ul[whichvar][PLUS0]   = U[whichvar][PLUS0]*ftilde + Ul[whichvar][PLUS0]*(1.0-ftilde);

      // Then monotonize
      monotonize(U[whichvar][PLUS0],Ur[whichvar][PLUS0],Ul[whichvar][PLUS0]);

      out_prims_r[whichvar].gf[index_arr[flux_dirn][PLUS0]] = Ur[whichvar][PLUS0];
      out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]] = Ul[whichvar][PLUS0];
    }
    // Ur depends on ftilde, which depends on points of U between MINUS2 and PLUS2
    out_prims_r[whichvar].gz_lo[flux_dirn]+=2;
    out_prims_r[whichvar].gz_hi[flux_dirn]+=2;
    // Ul depends on ftilde, which depends on points of U between MINUS2 and PLUS2
    out_prims_l[whichvar].gz_lo[flux_dirn]+=2;
    out_prims_l[whichvar].gz_hi[flux_dirn]+=2;
  }

  // *** LOOP 4: SHIFT Ur AND Ul ***
  /* Currently face values are set so that
   *      a) Ur(i) represents U(i+1/2), and
   *      b) Ul(i) represents U(i-1/2)
   *    Here, we shift so that the indices are consistent:
   *      a) U(i-1/2+epsilon) = oldUl(i)   = newUr(i)
   *      b) U(i-1/2-epsilon) = oldUr(i-1) = newUl(i)
   *    Note that this step is not strictly necessary if you keep
   *      track of indices when computing the flux. */
  for(int ww=0;ww<num_prims_to_reconstruct;ww++) {
    int whichvar=which_prims_to_reconstruct[ww];
    LOOP_DEFINE(3,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
      SET_INDEX_ARRAYS(-1,0,flux_dirn);
      temporary[index_arr[flux_dirn][PLUS0]] = out_prims_r[whichvar].gf[index_arr[flux_dirn][MINUS1]];
    }

    LOOP_DEFINE(3,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
      SET_INDEX_ARRAYS(0,0,flux_dirn);
      // Then shift so that Ur represents the gridpoint at i-1/2+epsilon,
      //                and Ul represents the gridpoint at i-1/2-epsilon.
      // Ur(i-1/2) = Ul(i-1/2)     = U(i-1/2+epsilon)
      // Ul(i-1/2) = Ur(i+1/2 - 1) = U(i-1/2-epsilon)
      out_prims_r[whichvar].gf[index_arr[flux_dirn][PLUS0]] = out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]];
      out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]] = temporary[index_arr[flux_dirn][PLUS0]];
    }
    // Ul was just shifted, so we lost another ghostzone.
    out_prims_l[whichvar].gz_lo[flux_dirn]+=1;
    out_prims_l[whichvar].gz_hi[flux_dirn]+=0;
    // As for Ur, we didn't need to get rid of another ghostzone,
    //    but we did ... seems wasteful!
    out_prims_r[whichvar].gz_lo[flux_dirn]+=1;
    out_prims_r[whichvar].gz_hi[flux_dirn]+=0;

  }
}

// Set SLOPE_LIMITER_COEFF = 2.0 for MC, 1 for minmod
#define SLOPE_LIMITER_COEFF 2.0

//Eq. 60 in JOURNAL OF COMPUTATIONAL PHYSICS 123, 1-14 (1996)
//   [note the factor of 2 missing in the |a_{j+1} - a_{j}| term].
//   Recall that dU = U_{i} - U_{i-1}.
static inline CCTK_REAL slope_limit(CCTK_REAL dU,CCTK_REAL dUp1) {
  if(dU*dUp1 > 0.0) {
    //delta_m_U=0.5 * [ (u_(i+1)-u_i) + (u_i-u_(i-1)) ] = (u_(i+1) - u_(i-1))/2  <-- first derivative, second-order; this should happen most of the time (smooth flows)
    CCTK_REAL delta_m_U = 0.5*(dU + dUp1);
    // EXPLANATION OF BELOW LINE OF CODE.
    // In short, sign_delta_a_j = sign(delta_m_U) = (0.0 < delta_m_U) - (delta_m_U < 0.0).
    //    If delta_m_U>0, then (0.0 < delta_m_U)==1, and (delta_m_U < 0.0)==0, so sign_delta_a_j=+1
    //    If delta_m_U<0, then (0.0 < delta_m_U)==0, and (delta_m_U < 0.0)==1, so sign_delta_a_j=-1
    //    If delta_m_U==0,then (0.0 < delta_m_U)==0, and (delta_m_U < 0.0)==0, so sign_delta_a_j=0
    int sign_delta_m_U = (0.0 < delta_m_U) - (delta_m_U < 0.0);
    //Decide whether to use 2nd order derivative or first-order derivative, limiting slope.
    return sign_delta_m_U*MIN(fabs(delta_m_U),MIN(SLOPE_LIMITER_COEFF*fabs(dUp1),SLOPE_LIMITER_COEFF*fabs(dU)));
  }
  return 0.0;
}

// standard Colella-Woodward parameters:
//    K0 = 0.1d0, eta1 = 20.0, eta2 = 0.05, epsilon = 0.01d0
#define K0      0.1
#define ETA1   20.0
#define ETA2    0.05
#define EPSILON 0.01
static inline void steepen_rho(const int eos_key,const int keytemp,const CCTK_REAL c2pprec,CCTK_REAL U[MAXNUMVARS][MAXNUMINDICES],CCTK_REAL slope_lim_dU[MAXNUMVARS][MAXNUMINDICES],CCTK_REAL gamma_therma,CCTK_REAL P_cold,CCTK_REAL gamma_cold,
                               CCTK_REAL *rho_br_ppm,CCTK_REAL *rho_bl_ppm) {

  DECLARE_CCTK_PARAMETERS;


  // Next compute centered differences d RHOB and d^2 RHOB
  CCTK_REAL d1rho_b     = 0.5*(U[RHOB][PLUS1] - U[RHOB][MINUS1]);
  CCTK_REAL d2rho_b_m1  = U[RHOB][PLUS0] - 2.0*U[RHOB][MINUS1] + U[RHOB][MINUS2];
  CCTK_REAL d2rho_b_p1  = U[RHOB][PLUS2] - 2.0*U[RHOB][PLUS1]  + U[RHOB][PLUS0];

  // Compute effective Gamma = (partial P / partial rho0)_s /(P/rho0)
  CCTK_REAL Gamma;
  if (!use_ConservativeToPrimitive){
    Gamma = gamma_therma + (gamma_cold-gamma_therma)*P_cold/U[PRESSURE][PLUS0];
  }else{
    CCTK_REAL xrho   = U[RHOB][PLUS0];
    CCTK_REAL xeps   = U[EPS][PLUS0];
    CCTK_REAL xye    = U[YE][PLUS0];
    CCTK_REAL xtemp  = U[TEMP][PLUS0];
    CCTK_REAL xpress = U[PRESSURE][PLUS0];
    CCTK_INT  keytemp = 0;     // We enter with eps
    CCTK_INT  npoints  = 1;
    CCTK_INT  anyerr, keyerr;

    EOS_Omni_Gamma(eos_key, keytemp, c2pprec, npoints, &xrho,&xtemp, &xye,  &xeps, &xpress, &Gamma, &keyerr, &anyerr);

    // printf("rhob %g, eps %g, ye %g, temp %g, press %g,  Gamma %g, Gamma poly %g err %d\n", xrho, xeps, xye, xtemp, xpress, Gamma,
    //        gamma_therma + (gamma_cold-gamma_therma)*P_cold/U[PRESSURE][PLUS0], anyerr);
  }


  CCTK_REAL contact_discontinuity_check = Gamma*K0*fabs(U[RHOB][PLUS1]-U[RHOB][MINUS1])*
    MIN(U[PRESSURE][PLUS1],U[PRESSURE][MINUS1])
    -fabs(U[PRESSURE][PLUS1]-U[PRESSURE][MINUS1])*MIN(U[RHOB][PLUS1],U[RHOB][MINUS1]);
  CCTK_REAL second_deriv_check = -d2rho_b_p1*d2rho_b_m1;
  CCTK_REAL relative_change_check = fabs(2.0*d1rho_b) - EPSILON*MIN(U[RHOB][PLUS1],U[RHOB][MINUS1]);

  if(contact_discontinuity_check >= 0.0 && second_deriv_check >= 0.0
     && relative_change_check >= 0.0) {

    CCTK_REAL eta_tilde=0.0;
    if (fabs(d1rho_b) > 0.0) {
      eta_tilde = -(1.0/6.0)*(d2rho_b_p1-d2rho_b_m1)/(2.0*d1rho_b);
    }
    CCTK_REAL eta = MAX(0.0,MIN(ETA1*(eta_tilde - ETA2),1.0));
    // Next compute Urp1 and Ul for RHOB, using the MC prescription:
    // Ur_p1 = U_p1   - 0.5*slope_lim_dU_p1
    CCTK_REAL rho_br_mc_p1 = U[RHOB][PLUS1] - 0.5*slope_lim_dU[RHOB][PLUS1];
    // Ul = U_m1 + 0.5*slope_lim_dU_m1
    // Based on this line of code, Ur[index] = a_j - \delta_m a_j / 2. (cf. Eq. 65 in Marti & Muller's "PPM Method for 1D Relativistic Hydro." paper)
    //    So: Ur[indexp1] = a_{j+1} - \delta_m a_{j+1} / 2. This is why we have rho_br_mc[indexp1]
    CCTK_REAL rho_bl_mc    = U[RHOB][MINUS1] + 0.5*slope_lim_dU[RHOB][MINUS1];

    rho_bl_ppm[PLUS0] = rho_bl_ppm[PLUS0]*(1.0-eta) + rho_bl_mc*eta;
    rho_br_ppm[PLUS0] = rho_br_ppm[PLUS0]*(1.0-eta) + rho_br_mc_p1*eta;

  }
}

static inline void monotonize(CCTK_REAL U,CCTK_REAL &Ur,CCTK_REAL &Ul) {
  CCTK_REAL dU = Ur - Ul;
  CCTK_REAL mU = 0.5*(Ur+Ul);

  if ( (Ur-U)*(U-Ul) <= 0.0) {
    Ur = U;
    Ul = U;
    return;
  }
  if ( dU*(U-mU) > (1.0/6.0)*SQR(dU)) {
    Ul = 3.0*U - 2.0*Ur;
    return;
  }
  if ( dU*(U-mU) < -(1.0/6.0)*SQR(dU)) {
    Ur = 3.0*U - 2.0*Ul;
    return;
  }
}

static inline void compute_P_cold__gamma_cold(CCTK_REAL rho_b,eos_struct &eos,   CCTK_REAL &P_cold,CCTK_REAL &gamma_cold) {
  // This code handles equations of state of the form defined
  // in Eqs 13-16 in http://arxiv.org/pdf/0802.0200.pdf

  // Default in case rho_b == 0.0
  if(rho_b==0.0) { P_cold = 0.0; gamma_cold = eos.gamma_tab[0]; return; }

  if(eos.neos==1) {
    // Eq. 14 of http://arxiv.org/pdf/0802.0200.pdf :
    // P_{cold} = K_i rho_i^{\Gamma_i}
    P_cold = eos.k_tab[0]*fasterpow_ppm_reconstruct(rho_b,eos.gamma_tab[0]);
    gamma_cold = eos.gamma_tab[0];
    return;
  }
  for(int nn=1;nn<eos.neos;nn++) {
    if (rho_b <= eos.rho_tab[nn]  &&  rho_b > eos.rho_tab[nn-1]) {
      P_cold = eos.k_tab[nn]*fasterpow_ppm_reconstruct(rho_b,eos.gamma_tab[nn]);
      gamma_cold = eos.gamma_tab[nn];
    }
  }
  if (rho_b > eos.rho_tab[eos.neos-1]) {
    P_cold = eos.k_tab[eos.neos]*fasterpow_ppm_reconstruct(rho_b,eos.gamma_tab[eos.neos]);
    gamma_cold = eos.gamma_tab[eos.neos];
  }
}

#define OMEGA1   0.75
#define OMEGA2  10.0
#define EPSILON2 0.33
static void ftilde_gf_compute(const cGH *cctkGH,const int *cctk_lsh,const int flux_dirn,gf_and_gz_struct *in_prims,CCTK_REAL *ftilde_gf) {
  int ijkgz_lo_hi[4][2];
  CCTK_REAL U[MAXNUMVARS][MAXNUMINDICES];
  /*Remove gcc unused variable warning/error Re: Pragma statement in loop define:*/
  CCTK_REAL dU,slope_lim_dU,Ur,Ul; dU=slope_lim_dU=Ur=Ul=0.0; dU*=0;
  // Compute ftilde, which is used for flattening left and right face values
  LOOP_DEFINE(2,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[VX+(flux_dirn-1)].gz_lo,in_prims[VX+(flux_dirn-1)].gz_hi) {
    SET_INDEX_ARRAYS(-2,2,flux_dirn);
    for(int ii=MINUS2;ii<=PLUS2;ii++) U[PRESSURE][ii] = in_prims[PRESSURE].gf[index_arr[flux_dirn][ii]];
    U[VX+(flux_dirn-1)][MINUS1] = in_prims[VX+(flux_dirn-1)].gf[index_arr[flux_dirn][MINUS1]];
    U[VX+(flux_dirn-1)][PLUS1]  = in_prims[VX+(flux_dirn-1)].gf[index_arr[flux_dirn][PLUS1]];

    // Compute ftilde, which is used for flattening left and right face values
    // DEPENDENCIES: P(MINUS2,MINUS1,PLUS1,PLUS2) and v^m(MINUS1,PLUS1), where m=flux_dirn={1,2,3}={x,y,z}.
    ftilde_gf[index_arr[flux_dirn][PLUS0]] = ftilde_compute(flux_dirn,U);
  }
}

static inline CCTK_REAL ftilde_compute(const int flux_dirn,CCTK_REAL U[MAXNUMVARS][MAXNUMINDICES]) {
  CCTK_REAL dP1 = U[PRESSURE][PLUS1] - U[PRESSURE][MINUS1];
  CCTK_REAL dP2 = U[PRESSURE][PLUS2] - U[PRESSURE][MINUS2];

  // MODIFICATION TO STANDARD PPM:
  // Cure roundoff error issues when dP1==0 or dP2==0 to 15 or more significant digits.
  CCTK_REAL avg1=0.5*(U[PRESSURE][PLUS1] + U[PRESSURE][MINUS1]);
  CCTK_REAL avg2=0.5*(U[PRESSURE][PLUS2] + U[PRESSURE][MINUS2]);
  if(fabs(dP1)/avg1<1e-15) dP1=0.0; /* If this is triggered, there is NO shock */
  if(fabs(dP2)/avg2<1e-15) dP2=0.0; /* If this is triggered alone, there may be a shock. Otherwise if triggered with above, NO shock. */

  CCTK_REAL dP1_over_dP2=1.0;
  if (dP2 != 0.0) dP1_over_dP2 = dP1/dP2;

  CCTK_REAL q1 = (dP1_over_dP2-OMEGA1)*OMEGA2;
  CCTK_REAL q2 = fabs(dP1)/MIN(U[PRESSURE][PLUS1],U[PRESSURE][MINUS1]);

  // w==0 -> NOT inside a shock
  CCTK_REAL w=0.0;

  // w==1 -> inside a shock
  if (q2 > EPSILON2 && q2*( (U[VX+(flux_dirn-1)][MINUS1]) - (U[VX+(flux_dirn-1)][PLUS1]) ) > 0.0) w = 1.0;

  return MIN(1.0, w*MAX(0.0,q1));
}
