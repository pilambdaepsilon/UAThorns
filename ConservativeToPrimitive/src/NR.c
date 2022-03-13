/*@@
   @file      NR.c
   @date      Sep 17, 2016
   @author    Daniel Siegel, Philipp Moesta
   @desc
   Routines for 3D Newton-Raphson root finding according
   to the scheme by Cerda-Duran (2008), Eq. 21, 22, 28 and 21, 22, 27.
   @enddesc
 @@*/


#include "con2prim.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

extern struct c2p_steer c2p;


void calc_prim_from_x(const double S_squared, const double BdotS,
      const double B_squared, const double * con, double * prim, double * x,
      const double g_cov[4][4]){

  /* Recover the primitive variables from the scalars (W,Z)
   * and conserved variables
   *
   * Eq (23)-(25) in Cerdá-Durán et al. 2008
   */

  double W = x[0];
  double Z = x[1];
  double T = x[2];

  double B1_cov, B2_cov, B3_cov;
  // Lower indices - covariant
  B1_cov = g_cov[1][1]*con[B1_con]+g_cov[1][2]*con[B2_con]+g_cov[1][3]*con[B3_con];
  B2_cov = g_cov[1][2]*con[B1_con]+g_cov[2][2]*con[B2_con]+g_cov[2][3]*con[B3_con];
  B3_cov = g_cov[1][3]*con[B1_con]+g_cov[2][3]*con[B2_con]+g_cov[3][3]*con[B3_con];

  prim[RHO] = con[D]/W;
  prim[v1_cov] = (con[S1_cov] + (BdotS)*B1_cov/Z)/(Z+B_squared);
  prim[v2_cov] = (con[S2_cov] + (BdotS)*B2_cov/Z)/(Z+B_squared);
  prim[v3_cov] = (con[S3_cov] + (BdotS)*B3_cov/Z)/(Z+B_squared);

  // calculate press, eps etc. from (rho, temp, Ye) using EOS
  // in this way, press and eps etc. are consistent with the other primitives
  const int keytemp=1;
  int keyerr=0;
  int anyerr=0;

  double xrho, xtemp, xye, xeps, xprs;
  xrho = prim[RHO];
  xtemp = T;
  xye = con[YE]/con[D];
  xeps = 0.0;
  xprs = 0.0;

  if (c2p.evolve_T){
    double xent = 0.0;
    double xabar = 0.0;
    EOS_press_ent_abar(c2p.eoskey,keytemp,xrho,&xeps,&xtemp,xye,&xprs,&xent,&xabar,&keyerr);
    prim[ENT] = xent;
    prim[A_BAR] = xabar;
    prim[MU_HAT] = 0.0;
  } else {
    EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &xrho, &xeps, &xtemp, &xye, &xprs, &keyerr, &anyerr);
  }

  prim[EPS] = xeps;
  prim[B1_con] = con[B1_con];
  prim[B2_con] = con[B2_con];
  prim[B3_con] = con[B3_con];
  prim[TEMP] = T;
  prim[YE] = xye;
  prim[PRESS] = xprs;
  prim[WLORENTZ] = W;

}


void calc_WZT_from_prim(double * prim, const double * con,
  const double g_con[4][4], double * x){

  // compute Newton-Raphson state vector x from primitive variables
  // (only needed for initial Newton-Raphson step)

  double W;

  /* compute W from 3-velocity components v^i (not recommended, using lorentz factor seems more robust)
  double v1_con = g_con[1][1]*prim[v1_cov]+g_con[1][2]*prim[v2_cov]+g_con[1][3]*prim[v3_cov];
  double v2_con = g_con[1][2]*prim[v1_cov]+g_con[2][2]*prim[v2_cov]+g_con[2][3]*prim[v3_cov];
  double v3_con = g_con[1][3]*prim[v1_cov]+g_con[2][3]*prim[v2_cov]+g_con[3][3]*prim[v3_cov];

  // v^2 = v^i * v_i
  double v_squared = prim[v1_cov]*v1_con + prim[v2_cov]*v2_con + prim[v3_cov]*v3_con;

  if(v_squared >= 1.0){
    W = 1.0e4;
  }else{
    W = 1.0/sqrt(1.0-v_squared);
  }
  */

  W = MAX(prim[WLORENTZ],1.0);

  // Always calculate rho from D and W
  // for consistency with other primitives
  prim[RHO] = con[D]/W;

  double P = P_from_prim(prim);
  double h = 1.0 + prim[EPS]+P/prim[RHO];
  double Z = prim[RHO]*h*W*W;

  x[0] = W;
  x[1] = Z;
  x[2] = prim[TEMP];
}


void calc_WZT_max(const double * con, const double B_squared, double * xmax,
                  int * have_xmax){

  // calculate maximum values for x = (W, Z, T) ("safe guess" initial values)
  // cf. Cerda-Duran et al. 2008, Eq. (39)-(42)

  double rhomax, epsmax, xtemp, xye, xprs;
  const int keytemp=0;
  int keyerr=0;
  int anyerr=0;

  rhomax = con[D];
  epsmax = (con[TAU] - B_squared/2.0) / con[D];

  // ensure that rhomax and epsmax are in validity range of EOS
  if (rhomax > c2p.eos_rho_max) {
    rhomax = 0.95*c2p.eos_rho_max;
  }
  if (epsmax > c2p.eos_eps_max) {
    epsmax = 0.95*c2p.eos_eps_max;
  }

  // compute Pmax, Tmax
  xtemp = 0.9*c2p.eos_temp_max; // initial guess, choose large enough
  xye = con[YE]/con[D];
  xprs = 0.0;

  EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &rhomax, &epsmax, &xtemp, &xye, &xprs, &keyerr, &anyerr);

  if(keyerr!=0){
    *have_xmax=0;
  }

  xmax[0] = 1.0e4;
  xmax[1] = con[TAU] + xprs + con[D] - 0.5*B_squared;
  xmax[2] = xtemp;

}


void NR_step_3D_P(double S_squared, double BdotS, double B_squared, const double * con,
  double * x, double dx[3], double f[3], int stepsize)
{
  /* Finding the roots of f(x):
   *
   * x_{n+1} = x_{n} - f(x)/J = x_{n} + dx_{n}
   *
   * where J is the Jacobian matrix J_{ij} = df_i/dx_j
   *
   * Here, compute dx = [dW, dZ, dT]
   */

  double J[3][3];
  double W = x[0];
  double Z = x[1];
  double T = x[2];

  double P,P_EOS,eps_EOS,dPdW,dPdZ,dPdT,dpEOSdrho,dpEOSdT;

  // compute partial derivatives of specific internal energy: dEdW and dEdZ
  EP_dPdW_dPdZ_dPdT(&eps_EOS,&P_EOS,&dPdW,&dPdZ,&dPdT,&dpEOSdrho,&dpEOSdT,x,con,stepsize);

  // compute pressure from state vector x and conservatives
  P = Z/(W*W) - con[D] * ( 1.0 + eps_EOS) / W ;

  // d/dW (1)
  J[0][0] = 2.0*((Z + B_squared)*(Z + B_squared)-S_squared-(2.0*Z + B_squared)*((BdotS)*(BdotS))/(Z*Z))*W;
  double a = J[0][0];

  // d/dZ (1)
  J[0][1] = (2.0*(Z + B_squared) + (2.0/(Z*Z) + 2.0*B_squared/(Z*Z*Z))*(BdotS*BdotS))*W*W - 2.0*(Z + B_squared);
  double b = J[0][1];

  // d/dT (1)
  J[0][2] = 0;
  double c = J[0][2];

  // d/dW (2)
  J[1][0] = 2.0*(con[TAU] + con[D] - Z - B_squared + (BdotS*BdotS)/(2.0*Z*Z) + P_EOS)*W - con[D]*dpEOSdrho;
  double d = J[1][0];

  // d/dZ (2)
  J[1][1] = (-1.0-(BdotS*BdotS)/(Z*Z*Z))*W*W;
  double e = J[1][1];

  // d/dT (2)
  J[1][2] = W*W*dpEOSdT;
  double fc = J[1][2];

  // d/dW (P-P(rho,T,Y_e))
  J[2][0] = dPdW;
  double g = J[2][0];

  // d/dZ (E-E(rho,T,Y_e))
  J[2][1] = dPdZ;
  double h = J[2][1];

  // d/dT (E-E(rho,T,Y_e))
  J[2][2] = dPdT;
  double k = J[2][2];

  // compute f(x) from (21), (22) and (28) in Cerda-Duran et al. 2008
  f[0] = ((Z+B_squared)*(Z+B_squared) - S_squared - (2.0*Z +B_squared)*(BdotS*BdotS)/(Z*Z))*W*W - ((Z + B_squared)*(Z + B_squared));
  f[1] = (con[TAU] + con[D] - Z - B_squared + (BdotS*BdotS)/(2.0*Z*Z) + P_EOS)*W*W + 0.5*B_squared;
  f[2] = P - P_EOS;

  // Compute the determinant
  double A = e*k-fc*h;
  double B = fc*g-d*k;
  double C = d*h-e*g;
  double detJ = a*(A) + b*(B) + c*(C);

  //Compute the matrix inverse
  double Ji[3][3];
  Ji[0][0] = A/detJ;
  Ji[1][0] = B/detJ;
  Ji[2][0] = C/detJ;
  Ji[0][1] = (c*h-b*k)/detJ;
  Ji[1][1] = (a*k-c*g)/detJ;
  Ji[2][1] = (g*b-a*h)/detJ;
  Ji[0][2] = (b*fc-c*e)/detJ;
  Ji[1][2] = (c*d-a*fc)/detJ;
  Ji[2][2] = (a*e-b*d)/detJ;

  // Compute the step size
  dx[0] = -Ji[0][0]* f[0] - Ji[0][1]* f[1] - Ji[0][2]* f[2];
  dx[1] = -Ji[1][0]* f[0] - Ji[1][1]* f[1] - Ji[1][2]* f[2];
  dx[2] = -Ji[2][0]* f[0] - Ji[2][1]* f[1] - Ji[2][2]* f[2];

}


void NR_step_3D_eps(double S_squared, double BdotS, double B_squared, const double * con,
  double * x, double dx[3], double f[3], int stepsize)
{
  /* Finding the roots of f(x):
   *
   * x_{n+1} = x_{n} - f(x)/J = x_{n} + dx_{n}
   *
   * where J is the Jacobian matrix J_{ij} = df_i/dx_j
   *
   * Here, compute dx = [dW, dZ, dT]
   */

  double J[3][3];
  double W = x[0];
  double Z = x[1];
  double T = x[2];

  double E,E_EOS,P_EOS,dEdW,dEdZ,dEdT,dpEOSdrho,dpEOSdT;

  // compute partial derivatives of specific internal energy: dEdW and dEdZ
  EP_dEdW_dEdZ_dEdT(&E_EOS,&P_EOS,&dEdW,&dEdZ,&dEdT,&dpEOSdrho,&dpEOSdT,x,con,stepsize);
  printf("E: %g, %g, %g, %g, %g, %g, %g, %g\n", E,E_EOS,P_EOS,dEdW,dEdZ,dEdT,dpEOSdrho,dpEOSdT);

  // compute specific internal energy from state vector x and conservatives
  // Eq. (25) in Cerda-Duran et al. 2008
  E = -1.0 + Z/(con[D]*W) - P_EOS*W/con[D];
  E=fmax(E,c2p.eos_eps_min);

  // d/dW (1)
  J[0][0] = 2.0*((Z + B_squared)*(Z + B_squared)-S_squared-(2.0*Z + B_squared)*((BdotS)*(BdotS))/(Z*Z))*W;
  double a = J[0][0];

  // d/dZ (1)
  J[0][1] = (2.0*(Z + B_squared) + (2.0/(Z*Z) + 2.0*B_squared/(Z*Z*Z))*(BdotS*BdotS))*W*W - 2.0*(Z + B_squared);
  double b = J[0][1];

  // d/dT (1)
  J[0][2] = 0;
  double c = J[0][2];

  // d/dW (2)
  J[1][0] = 2.0*(con[TAU] + con[D] - Z - B_squared + (BdotS*BdotS)/(2.0*Z*Z) + P_EOS)*W - con[D]*dpEOSdrho;
  double d = J[1][0];

  // d/dZ (2)
  J[1][1] = (-1.0-(BdotS*BdotS)/(Z*Z*Z))*W*W;
  double e = J[1][1];

  // d/dT (2)
  J[1][2] = W*W*dpEOSdT;
  double fc = J[1][2];

  // d/dW (E-E(rho,T,Y_e))
  J[2][0] = dEdW;
  double g = J[2][0];

  // d/dZ (E-E(rho,T,Y_e))
  J[2][1] = dEdZ;
  double h = J[2][1];

  // d/dT (E-E(rho,T,Y_e))
  J[2][2] = dEdT;
  double k = J[2][2];

  // compute f(x) from (21), (22) and (28) in Cerda-Duran et al. 2008
  f[0] = ((Z+B_squared)*(Z+B_squared) - S_squared - (2.0*Z +B_squared)*(BdotS*BdotS)/(Z*Z))*W*W - ((Z + B_squared)*(Z + B_squared));
  f[1] = (con[TAU] + con[D] - Z - B_squared + (BdotS*BdotS)/(2.0*Z*Z) + P_EOS)*W*W + 0.5*B_squared;
  f[2] = E - E_EOS;

  // Compute the determinant
  double A = e*k-fc*h;
  double B = fc*g-d*k;
  double C = d*h-e*g;
  double detJ = a*(A) + b*(B) + c*(C);

  //Compute the matrix inverse
  double Ji[3][3];
  Ji[0][0] = A/detJ;
  Ji[1][0] = B/detJ;
  Ji[2][0] = C/detJ;
  Ji[0][1] = (c*h-b*k)/detJ;
  Ji[1][1] = (a*k-c*g)/detJ;
  Ji[2][1] = (g*b-a*h)/detJ;
  Ji[0][2] = (b*fc-c*e)/detJ;
  Ji[1][2] = (c*d-a*fc)/detJ;
  Ji[2][2] = (a*e-b*d)/detJ;

  // Compute the step size
  dx[0] = -Ji[0][0]* f[0] - Ji[0][1]* f[1] - Ji[0][2]* f[2];
  dx[1] = -Ji[1][0]* f[0] - Ji[1][1]* f[1] - Ji[1][2]* f[2];
  dx[2] = -Ji[2][0]* f[0] - Ji[2][1]* f[1] - Ji[2][2]* f[2];

}


void NR_3D(struct c2p_report * c2p_rep, const double S_squared, const double BdotS,
      const double B_squared, const double * con, double * prim, const double g_con[4][4],
      const double g_cov[4][4], const double tol_x, int SAFEGUESS, int stepsize){

  // 3D Newton-Raphson scheme, using state vector x = (W, Z, T) and 3D function
  // f(x) = (f1(x), f2(x), f3(x)) given by Eqs. (21), (22), (29) of Cerda-Duran et al. 2008

  // initialize Newton-Raphson state vector x  = (W, Z, T)
  double x[3];
  for (int l=0;l<3;l++){
    x[l] = 0.0;
  }

  double f[3];      // root finding function f
  double dx[3];     // displacement vector
  double x_lowlim[3]; // lower limits on W, Z, T
  double x_old[3];  // old state vector
  double error[3];  // error vector abs(dx/x)
  int done = 0;     // flag to control iteration loop
  int count = 0;    // count number of iterations

  // define lower limits for x
  x_lowlim[0] = 1.0;
  x_lowlim[1] = 0.0;
  x_lowlim[2] = c2p.eos_temp_min;

  //termination critera definitions
  double MIN_NEWT_TOL = tol_x;
  double NEWT_TOL=tol_x;
  int MAX_NEWT_ITER = c2p.max_iterations;
  int EXTRA_NEWT_ITER = c2p.extra_iterations;


  // compute initial guess for Newton-Raphson state vector x
  if(SAFEGUESS==0){
    calc_WZT_from_prim(prim,con,g_con,x);
  }
  else if(SAFEGUESS==1){
    int have_safe_guess=1;
    calc_WZT_max(con,B_squared,x,&have_safe_guess);
    if(have_safe_guess==0){
      printf("Cannot find safe initial guess! Giving up...\n");
      c2p_rep->failed = true;
      c2p_rep->count = count;
      return;
    }
  }

  // initialize variables
  for(int i=0;i<3;i++){
    x_old[i] = 0.0;
    dx[i] = 0.0;
    f[i] = 0.0;
    error[i] = 0.0;
    //check for NaNs
    if (x[i] != x[i]){
      printf("%e\t%e\t%e\n", x[0],x[1],x[2]);
      printf("3D_NR Error: NaN in initial NR guesses...\n");
      c2p_rep->failed = true;
      c2p_rep->count = count;
      return;
    }
  }


  // start NR scheme
  int i_extra;
  int doing_extra;
  i_extra = doing_extra = 0;

  if (EXTRA_NEWT_ITER==0) {
    i_extra=-1;
  }

  bool keep_iterating=1;
  double maxerror;

  while (keep_iterating) {

    // do Newton-Raphson step
    if (c2p.c2p_method == 1){
      NR_step_3D_eps(S_squared, BdotS, B_squared, con, x, dx, f, stepsize);
    } else if (c2p.c2p_method == 2) {
      NR_step_3D_P(S_squared, BdotS, B_squared, con, x, dx, f, stepsize);
    }

    printf("S_squared: %g, BdotS: %g, B_squared: %g\n", S_squared, BdotS, B_squared);


    // update x and exclude unphysical regimes
    for(int i=0;i<3;i++){
      x_old[i] = x[i];
      x[i] = fabs(x[i] + dx[i] - x_lowlim[i]) + x_lowlim[i];
      error[i] = fabs((x[i]-x_old[i])/x[i]);

      if (x[i] != x[i]){
        printf("3D_NR Error: NaN in NR state vector...\n");
        printf("%e\t%e\t%e\n", x[0],x[1],x[2]);
        printf("Iteration No. %d\n", count);
        printf("x_old[i]\tx_new[i]\t\tdx[i]\t\tf[i]\t\terror[i]\n");
        printf("%e\t%e\t%e\t%e\t%e\n", x_old[0],x[0],dx[0],f[0],error[0]);
        printf("%e\t%e\t%e\t%e\t%e\n", x_old[1],x[1],dx[1],f[1],error[1]);
        printf("%e\t%e\t%e\t%e\t%e\n", x_old[2],x[2],dx[2],f[2],error[2]);
        c2p_rep->failed = true;
        c2p_rep->count = count;
        return;
      }
    }
    maxerror = MAX(error[0], MAX(error[1],error[2]));

#if DEBUG
    // report current NR step
    printf("-----\n");
    printf("Iteration No. %d\n", count);
    printf("x_old[i]\tx_new[i]\t\tdx[i]\t\tf[i]\t\terror[i]\n");
    printf("%e\t%e\t%e\t%e\t%e\n", x_old[0],x[0],dx[0],f[0],error[0]);
    printf("%e\t%e\t%e\t%e\t%e\n", x_old[1],x[1],dx[1],f[1],error[1]);
    printf("%e\t%e\t%e\t%e\t%e\n", x_old[2],x[2],dx[2],f[2],error[2]);
    printf("maxerror = %15e\n", maxerror);
#endif

    ++count;

    // termination criterion

    if( (fabs(maxerror) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;


    if( ((fabs(maxerror) <= NEWT_TOL)&&(doing_extra == 0))
        || (i_extra >= EXTRA_NEWT_ITER) || (count >= (MAX_NEWT_ITER)) ) {
      keep_iterating = 0;
    }
  }

  //  Check for bad untrapped divergences
  if( (!isfinite(f[0])) || (!isfinite(f[1]))|| (!isfinite(f[2])) ) {
    c2p_rep->failed = true;
  } else if( fabs(maxerror) <= NEWT_TOL ){
    c2p_rep->failed = false;
  } else if( (fabs(maxerror) <= MIN_NEWT_TOL) && (fabs(maxerror) > NEWT_TOL) ){
    c2p_rep->failed = false;
  } else {
    c2p_rep->failed = true;
  }

  // Recover the primitive variables from the final scalars x = (W,Z,T)
  calc_prim_from_x(S_squared, BdotS, B_squared, con, prim, x, g_cov);
  c2p_rep->nEOScalls += 1;
  c2p_rep->count = count;
}
