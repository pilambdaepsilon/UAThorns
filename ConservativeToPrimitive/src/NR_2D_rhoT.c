#include "con2prim.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


void calc_rhoT_max(const double * con, const double B_squared, double * xmax,
                  int * have_xmax){

  // calculate maximum values for x = (rho, T) ("safe guess" initial values)
  // cf. Cerda-Duran et al. 2008, Eq. (39)-(42)

  double rhomax = con[D];
  double epsmax = (con[TAU] - B_squared/2.0) / con[D];

  // ensure that rhomax and epsmax are in validity range of EOS
  if (rhomax > c2p.eos_rho_max) {
    rhomax = 0.95*c2p.eos_rho_max;
  }
  if (epsmax > c2p.eos_eps_max) {
    epsmax = 0.95*c2p.eos_eps_max;
  }

  // compute Pmax, Tmax
  const int keytemp=0;
  const double prec = 1.0e-10;
  int keyerr=0;
  int anyerr=0;

  double xtemp, xye, xprs;
  xtemp = 0.9*c2p.eos_temp_max; // initial guess, choose large enough
  xye = con[YE]/con[D];
  xprs = 0.0;

  EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &rhomax, &epsmax, &xtemp, &xye, &xprs, &keyerr, &anyerr);

  if(keyerr!=0){
    *have_xmax=0;
  }

  xmax[0] = rhomax;
  xmax[1] = xtemp;

}


void calc_prim_from_x_2D_rhoT(const double S_squared, const double BdotS,
      const double B_squared, const double * con, double * prim, double * x,
      const double glo[4][4]){

  /* Recover the primitive variables from the scalars (W,Z)
   * and conserved variables, Eq. (23)-(25) in Cerdá-Durán et al. 2008
   */

  double rho = x[0];
  double T   = x[1];
  double W = con[D] / rho;

  // calculate press, eps etc. from (rho, temp, Ye) using EOS,
  // required for consistency
  const int keytemp=1;
  int keyerr=0;
  int anyerr=0;

  double xrho, xtemp, xye, xeps, xprs;
  xrho = rho;
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

  double Z = rho * (1.0 + xeps + xprs/rho) * W*W;

  double B1_cov, B2_cov, B3_cov;
  // Lower indices - covariant
  B1_cov = glo[1][1]*con[B1_con]+glo[1][2]*con[B2_con]+glo[1][3]*con[B3_con];
  B2_cov = glo[1][2]*con[B1_con]+glo[2][2]*con[B2_con]+glo[2][3]*con[B3_con];
  B3_cov = glo[1][3]*con[B1_con]+glo[2][3]*con[B2_con]+glo[3][3]*con[B3_con];

  prim[RHO] = rho;
  prim[v1_cov] = (con[S1_cov] + (BdotS)*B1_cov/Z)/(Z+B_squared);
  prim[v2_cov] = (con[S2_cov] + (BdotS)*B2_cov/Z)/(Z+B_squared);
  prim[v3_cov] = (con[S3_cov] + (BdotS)*B3_cov/Z)/(Z+B_squared);

  prim[EPS] = xeps;
  prim[B1_con] = con[B1_con];
  prim[B2_con] = con[B2_con];
  prim[B3_con] = con[B3_con];
  prim[TEMP] = T;
  prim[YE] = xye;
  prim[PRESS] = xprs;
  prim[WLORENTZ] = W;

}


void NR_step_2D_rhoT(double S_squared, double BdotS, double B_squared, const double * con,
		double * x, double dx[2], double f[2], int stepsize)
{
  /* Finding the roots of f(x):
   *
   * x_{n+1} = x_{n} - f(x)/J = x_{n} + dx_{n}
   *
   * where J is the Jacobian matrix J_{ij} = df_i/dx_j
   *
   * Here, compute dx = [drho, dT]
   */

  double J[2][2];
  double Ji[2][2];
  double rho = x[0];
  double T   = x[1];

  double P,E,H,W,Z,df1_dZ,df2_dZ,df1_dW,df2_dW, df2_dP, dZ_drho, dZ_dT, dW_drho;
  double dEdrho = 0.0;
  double dEdT = 0.0;
  double dPdrho = 0.0;
  double dPdT = 0.0;

  W  = con[D]/rho;

  /* need partial derivatives of specific internal energy and pressure wrt density and
   * temperature. Those need to be based on primitives computed from Newton-Raphson state
   * vector x and conservatives
   */
  EOS_EP_dEdr_dEdt_dPdr_dPdt_2D(rho,T,con,&E,&P,&dEdrho,&dEdT,&dPdrho,&dPdT,stepsize);

  H = 1.0 + E + P/rho;
  Z = rho * H * W*W;

  dW_drho = - W / rho;
  dZ_drho = - W*W*H + con[D]*W * ( dEdrho - P/(rho*rho) + dPdrho/rho );
  dZ_dT = con[D]*W * ( dEdT + dPdT/rho );

  df1_dZ = (2.0*(Z + B_squared) + (2.0/(Z*Z) + 2.0*B_squared/(Z*Z*Z))*(BdotS*BdotS))*W*W - 2.0*(Z + B_squared);
  df1_dW = 2.0*((Z + B_squared)*(Z + B_squared)-S_squared-(2.0*Z + B_squared)*((BdotS)*(BdotS))/(Z*Z))*W;

  df2_dZ = (-1.0-(BdotS*BdotS)/(Z*Z*Z))*W*W;
  df2_dW = 2.0*(con[TAU] + con[D] - Z - B_squared + (BdotS*BdotS)/(2.0*Z*Z) + P)*W;
  df2_dP = W*W;

  // df1_drho
  J[0][0] = df1_dZ * dZ_drho + df1_dW * dW_drho;
  double a = J[0][0];

  // df1_dT
  J[0][1] = df1_dZ * dZ_dT;
  double b = J[0][1];

  // df2_drho
  J[1][0] = df2_dZ * dZ_drho + df2_dW * dW_drho;
  double c = J[1][0];

  // df2_dT
  J[1][1] = df2_dZ * dZ_dT + df2_dP * dPdT;
  double d = J[1][1];


  // compute f(x) from (27), (28) in Siegel et al. 2018
  f[0] = ((Z+B_squared)*(Z+B_squared) - S_squared - (2.0*Z +B_squared)*(BdotS*BdotS)/(Z*Z))*W*W - ((Z + B_squared)*(Z + B_squared));
  f[1] = (con[TAU] + con[D] - Z - B_squared + (BdotS*BdotS)/(2.0*Z*Z) + P)*W*W + 0.5*B_squared;

  double detJ = a*d - b*c;
  Ji[0][0] = d/detJ;
  Ji[0][1] = -b/detJ;
  Ji[1][0] = -c/detJ;
  Ji[1][1] = a/detJ;

  // Compute the step size
  dx[0] = -Ji[0][0]* f[0] - Ji[0][1]* f[1];
  dx[1] = -Ji[1][0]* f[0] - Ji[1][1]* f[1];

}


void NR_2D_rhoT(struct c2p_report * c2p_rep, const double S_squared,
      const double BdotS, const double B_squared, const double * con, double * prim,
	  const double glo[4][4], const double tol_x, int SAFEGUESS, int stepsize){

  // 2D Newton-Raphson scheme, using state vector x = (rho, T) and 2D function
  // f(x) = (f1(x), f2(x)) given by Eqs. (27), (28) of Siegel et al. 2018

  // initialize Newton-Raphson state vector x  = (rho, T)
  double x[2];
  for (int l=0;l<2;l++){
    x[l] = 0.0;
  }
  double f[2];      // root finding function f
  double x_lowlim[2]; // lower limits on rho, T
  double dx[2];     // displacement vector
  double x_old[2];  // old state vector
  double error[2];  // error vector abs(dx/x)
  int done = 0;     // flag to control iteration loop
  int count = 0;    // count number of iterations

  x_lowlim[0] = c2p.eos_rho_min;
  x_lowlim[1] = c2p.eos_temp_min;

  //termination critera definitions
  double MIN_NEWT_TOL = tol_x;
  double NEWT_TOL=tol_x;
  int MAX_NEWT_ITER = c2p.max_iterations;
  int EXTRA_NEWT_ITER = c2p.extra_iterations;

  // set initial guess for Newton-Raphson state vector x = (rho, T)
  if(SAFEGUESS==0){
    x[0] = prim[RHO];
    x[1] = prim[TEMP];
  }
  else if(SAFEGUESS==1){
    int have_safe_guess=1;
    calc_rhoT_max(con,B_squared,x,&have_safe_guess);
    if(have_safe_guess==0){
#if DEBUG
      printf("Cannot find safe initial guess! Giving up...\n");
#endif
      c2p_rep->failed = true;
      c2p_rep->count = 0;
    }
  }

  // initialize variables
  int i;
  for(i=0;i<2;i++){
    x_old[i] = 0.0;
    dx[i] = 0.0;
    f[i] = 0.0;
    error[i] = 0.0;
    //check for NaNs
    if (x[i] != x[i]){
#if DEBUG2
      printf("%e\t%e\n", x[0],x[1]);
      printf("2D_NR Error: NaN in initial NR guesses...\n");
#endif
      c2p_rep->failed = true;
      c2p_rep->count = 0;
      return;
    }
  }

  int i_extra = 0;
  int doing_extra = 0;

  if (EXTRA_NEWT_ITER==0) {
    i_extra=-1;
  }

  bool keep_iterating=1;
  double maxerror;

  while (keep_iterating) {

    // do Newton-Raphson step
    NR_step_2D_rhoT(S_squared, BdotS, B_squared, con, x, dx, f, stepsize);

    // Update x vector and compute error
    for(i=0;i<2;i++){
      x_old[i] = x[i];
      // update x and exclude unphysical regime
      x[i] = fmax(x[i] + dx[i], x_lowlim[i]);
      error[i] = fabs((x[i]-x_old[i])/x[i]);

      if (x[i] != x[i]){
#if DEBUG2
        printf("2D_NR Error: NaN in initial NR guesses...\n");
        printf("%e\t%e\n", x[0],x[1]);
        printf("Iteration No. %d\n", count);
        printf("x_old[i]\tx_new[i]\t\tdx[i]\t\tf[i]\t\terror[i]\n");
        printf("%e\t%e\t%e\t%e\t%e\n", x_old[0],x[0],dx[0],f[0],error[0]);
        printf("%e\t%e\t%e\t%e\t%e\n", x_old[1],x[1],dx[1],f[1],error[1]);
        printf("maxerror = %15e\n", maxerror);
#endif
        c2p_rep->failed = true;
        c2p_rep->count = 0;
        return;
      }
    }
    maxerror = MAX(error[0], error[1]);

#if DEBUG
    // report current NR step
    printf("-----\n");
    printf("Iteration No. %d\n", count);
    printf("x_old[i]\tx_new[i]\t\tdx[i]\t\tf[i]\t\terror[i]\n");
    printf("%e\t%e\t%e\t%e\t%e\n", x_old[0],x[0],dx[0],f[0],error[0]);
    printf("%e\t%e\t%e\t%e\t%e\n", x_old[1],x[1],dx[1],f[1],error[1]);
    printf("maxerror = %15e\n", maxerror);
#endif

    count ++;

    // termination criterion
    if( (fabs(maxerror) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(maxerror) <= NEWT_TOL)&&(doing_extra == 0))
        || (i_extra >= EXTRA_NEWT_ITER) || (count >= (MAX_NEWT_ITER)) ) {
      keep_iterating = 0;
    }

  }  // END of while(keep_iterating)

  //  Check for bad untrapped divergences
  if( (!isfinite(f[0])) || (!isfinite(f[1])) ) {
    c2p_rep->failed = true;
  }

  if( fabs(maxerror) <= NEWT_TOL ){
    c2p_rep->failed = false;
  }
  else if( (fabs(maxerror) <= MIN_NEWT_TOL) && (fabs(maxerror) > NEWT_TOL) ){
    c2p_rep->failed = false;
  }
  else {
    c2p_rep->failed = true;
  }

  // Recover the primitive variables from the final scalars x = (rho,T)
  calc_prim_from_x_2D_rhoT(S_squared, BdotS, B_squared, con, prim, x, glo);
  c2p_rep->count = count;
}
