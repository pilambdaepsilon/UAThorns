#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include "con2prim.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#define W_TOO_BIG       (1.e20)    /* \gamma^2 (\rho_0 + u + p) is assumed
                                      to always be smaller than this.  This
                                      is used to detect solver failures */
#define SQR(x) ((x) * (x))
extern struct c2p_steer c2p;

struct LocGlob {
  double Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Dcon,half_Bsq ;
} ;

// function declarations:
//static double vsq_calc(double W, struct LocGlob *lgp);
//static double x1_of_x0(double x0, struct LocGlob *lgp ) ;

static int twod_newton_raphson( double x[], struct LocGlob *lgp, const double * con,
      double * temp_guess,
      void (*funcd) (double [], double [], double [], double [][2], double *, double *,
                     struct LocGlob *, const double *, double *,  const double, double),
      struct c2p_report * c2p_rep , const double tol_x);

static void func_vsq( double [], double [], double [], double [][2], double *f, double *df,
      struct LocGlob *lgp, const double * con, double * temp_guess,
      const double tol_x, double h_min);


static double eos_quantities_general2(const double * x, double *dpdw, double *dpdvsq,
      const struct LocGlob *lgp, const double * con, double * temp_guess,
      const double tol_x);


void comp_press_eps_from_x2(double * press, double * eps, const double * x, const double BdotS,
      const double B_squared, const double * con );

void calc_prim_harm(struct LocGlob *lgp, const double * con, double * prim, const double temp_guess, 
     double * x, const double g_cov[4][4], const double g_con[4][4], struct c2p_report * c2p_rep, const double * METRIC, const double * METRIC_PHYS, const double * METRIC_LAP_PSI4, 
     const int iteration, const double xco, const double yco, const double zco);


/////////////////////////////////////////////////////////////////

/******************************************************************************
    Main Con2Prim routine: NR_2D_Noble
 ******************************************************************************

  -- Inversion from GRMHD conserved variables to primitive variables assuming
     a guess based on given primitives in prim

  -- Uses the 2D method of Noble et al. 2006, ApJ 641, 626:
       -- Solves for two independent variables (W,v^2) via a 2D
          Newton-Raphson method
       -- This implementation is largely based on the publicly available implementation
          by Noble et al. 2006
       -- Can be used with a general 3-parameter EOS in this implementation

  -- Note that the notation used herein is largely that of Noble et al. (2006)


INPUT: (using notation in Siegel et al. 2018, ApJS)

  S_squared = S_i S^i = S^2
  BdotS = B_i * S^i
  B_squared = B_i * B^i
  con = array of conservatives
  prim = array of primitives (guesses)
  g_con = g^munu # contravariant components of the space-time metric
  g_cov = g_munu # covariant components of the space-time metric

OUTPUT:  (using GRHydro variable defintions)
  prim  = array of recovered primitives
  c2p_rep = con2prim structure containing variables for reporting back to con2prim framework
  Note: c2_rep.c2p_keyerr = (i*100 + j)  where
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used)
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence
                   (occurrence of "nan" or "+/-inf" ;

         j = 0 -> success
             1 -> failure: some sort of failure in Newton-Raphson;
             2 -> failure: unphysical vsq = v^2  value at initial guess;
             3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1
             ( used to be  5 -> failure: rho,uu <= 0   but now sets epsnegative to non-zero )

**********************************************************************************/
void HARM(struct c2p_report * c2p_rep, const double S_squared, const double BdotS,
      const double B_squared, const double * con, double * prim, CCTK_REAL g_con[4][4],
      CCTK_REAL g_cov[4][4], const double tol_x,double xco, double yco, double zco,const int iteration,const int index, double * METRIC, double * METRIC_PHYS, double * METRIC_LAP_PSI4)
{
  double x_2d[2];
  double sx, sy, sz;
  double usx, usy, usz;
  double tau, dens;
  double QdotB;
  double rho0,u,p,w,gammasq,gamma,gtmp,W_last,W,vsq;
  double g_o_WBsq, QdB_o_W;
  int  i_increase;

  struct LocGlob lg;

  c2p_rep->c2p_keyerr = 0;

  tau  = con[TAU];
  dens = con[D];

/*============================ START NEW ============================*/
  //From Illinois, if guess is already bad, then stop.
  if( con[D] <= 0. ){c2p_rep->failed=true; c2p_rep->c2p_keyerr=-100; return;}
/*============================ END NEW ============================*/

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  // NEW: All B-related quantities now carry appropriate 
  // factors of 1/sqrt(4*pi) to make consistent with Illinois
  lg.Bsq = B_squared;

  lg.QdotB = BdotS;
  lg.QdotBsq = lg.QdotB*lg.QdotB ;

  lg.Qdotn = -(tau + dens);

  lg.Qtsq = S_squared;

  lg.Dcon = dens;

  lg.half_Bsq = 0.5*lg.Bsq;

  /* calculate W from last timestep and use for guess -- changed this to use 3-metric */
  double velupx = METRIC_PHYS[C2P_GUPXX]*prim[v1_cov] + METRIC_PHYS[C2P_GUPXY]*prim[v2_cov] + METRIC_PHYS[C2P_GUPXZ]*prim[v3_cov];
  double velupy = METRIC_PHYS[C2P_GUPXY]*prim[v1_cov] + METRIC_PHYS[C2P_GUPYY]*prim[v2_cov] + METRIC_PHYS[C2P_GUPYZ]*prim[v3_cov];
  double velupz = METRIC_PHYS[C2P_GUPXZ]*prim[v1_cov] + METRIC_PHYS[C2P_GUPYZ]*prim[v2_cov] + METRIC_PHYS[C2P_GUPZZ]*prim[v3_cov];
  //double velupx = g_con[1][1]*prim[v1_cov] + g_con[1][2]*prim[v2_cov] + g_con[1][3]*prim[v3_cov];
  //double velupy = g_con[2][1]*prim[v1_cov] + g_con[2][2]*prim[v2_cov] + g_con[2][3]*prim[v3_cov];
  //double velupz = g_con[3][1]*prim[v1_cov] + g_con[3][2]*prim[v2_cov] + g_con[3][3]*prim[v3_cov];

  vsq = velupx*prim[v1_cov] + velupy*prim[v2_cov] + velupz*prim[v3_cov];

  if( (vsq < 0.) && (fabs(vsq) < 1.0e-13) ) {
    vsq = fabs(vsq);
  }
  if(vsq < 0. || vsq >= 1. ) {
    c2p_rep->c2p_keyerr = 2;
    return;
  }

  gammasq = 1. / (1. - vsq);
  gamma  = sqrt(gammasq);

  // always calculate rho from D and gamma, required for
  // consistency with other primitives
  rho0 = lg.Dcon / gamma ;
  double eps_last = fmax(prim[EPS], c2p.eos_eps_min);
  u = eps_last * rho0;

  const int keytemp = 1; // don't know the temperature, only eps
  int keyerr = 0;
  int anyerr = 0;

  double xeps = eps_last;
  double xtemp = prim[TEMP];
  double ye = con[YE]/con[D];
  double xP = 0.;

  EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &rho0, &xeps, &xtemp, &ye, &xP, &keyerr, &anyerr);

  p = xP;
  w = rho0 + u + p ;

  W_last = w*gammasq ; // this is z = rho h W^2 in NR 2D, 3D

  // Make sure that W is large enough so that v^2 < 1 :
  //i_increase = 0;
  //while( (( W_last*W_last*W_last * ( W_last + 2.*lg.Bsq )
  //          - lg.QdotBsq*(2.*W_last + lg.Bsq) ) <= W_last*W_last*(lg.Qtsq-lg.Bsq*lg.Bsq))
  //       && (i_increase < 10) ) {
  //  W_last *= 10.;
  //  i_increase++;
  //}
  // Calculate W and vsq:
  //x_2d[0] = fabs( W_last );
  //x_2d[1] = x1_of_x0( W_last, &lg ) ;


  // initialize W and vsq:
  // note: v^2 < 1 is already ensured
  x_2d[0] = fabs( W_last );
  x_2d[1] = vsq;

  // temperature guess
  double temp_guess = prim[TEMP];

#if DEBUG
  printf("*****************************************************************\n");
  printf("initial rho = %e\n", rho0);
  printf("initial eps = %e\n", eps_last);
  printf("initial press = %e\n", p);
  printf("initial enth = %e\n", 1.0+eps_last + p/rho0);
  printf("initial W=rho*h*W^2 = %e\n", W_last);
  printf("initial vsq = %e\n", vsq);
  printf("initial gamma = %e\n", gamma);
  printf("x_2d initial= %e, %e\n", x_2d[0], x_2d[1]);
  printf("temp_guess= %e\n", temp_guess);
#endif

  // start recovery
  c2p_rep->c2p_keyerr = twod_newton_raphson(x_2d, &lg, con, &temp_guess, func_vsq, c2p_rep, tol_x) ;

  W = x_2d[0];
  vsq = x_2d[1];

  /* Return in case of failure */
  if( c2p_rep->c2p_keyerr != 0.) {
    c2p_rep->c2p_keyerr = c2p_rep->c2p_keyerr*100.+1.;
    return;
  }
  else{
    if(W <= 0. || W > W_TOO_BIG) {
      c2p_rep->c2p_keyerr = c2p_rep->c2p_keyerr*100.+3.;
      return;
    }
  }

  if( vsq >= 1. ) {
    c2p_rep->c2p_keyerr = c2p_rep->c2p_keyerr*100.+4.;
    return;
  }

  // compute final set of primitives
  calc_prim_harm(&lg, con, prim, temp_guess, x_2d, g_cov, g_con, c2p_rep, METRIC, METRIC_PHYS, METRIC_LAP_PSI4, iteration, xco, yco, zco);

  return;
}

/******************************************************************************
    Further routines:
 ******************************************************************************/


/****************************************************************************
   vsq_calc():

      -- evaluate v^2 (spatial, normalized velocity) from
            W = \gamma^2 w

****************************************************************************/
/*
static double vsq_calc(double W, struct LocGlob *lgp)
{
        double Wsq,Xsq,Bsq_W;

        Wsq = W*W ;
        Bsq_W = (lgp->Bsq + W);
        Xsq = Bsq_W * Bsq_W;

        return(  ( Wsq * lgp->Qtsq  + lgp->QdotBsq * (Bsq_W + W)) / (Wsq*Xsq) );
}
*/

/********************************************************************

  x1_of_x0():

    -- calculates v^2 from W  with some physical bounds checking;
    -- asumes x0 is already physical
    -- makes v^2 physical  if not;

*********************************************************************/
/*
static double x1_of_x0(double x0, struct LocGlob *lgp )
{
  double vsq;
  double dv = 1.e-15;

  vsq = fabs(vsq_calc(x0,lgp)) ; // guaranteed to be positive

  return( ( vsq >= 1. ) ? (1.0 - dv) : vsq   );

}
*/

/************************************************************

  twod_newton_raphson():

    -- performs Newton-Rapshon method on an 2D system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static int twod_newton_raphson( double x[], struct LocGlob *lgp, const double * con,
      double * temp_guess,
      void (*funcd) (double [], double [], double [], double [][2], double *, double *,
                     struct LocGlob *, const double *, double *, const double, double),
      struct c2p_report * c2p_rep, const double tol_x)
{
  double f, df, dx[2], x_old[2];
  double resid[2], jac[2][2];
  double errx, x_orig[2];
  int    n_iter, i_extra, doing_extra;
  double vsq,W,W_old;
  const double dv = (1.-1.e-15);
  dx[0]=dx[1]=0.0;

  double NEWT_TOL = tol_x;
  double MIN_NEWT_TOL = tol_x;
  double MAX_NEWT_ITER = c2p.max_iterations;
  double EXTRA_NEWT_ITER = c2p.extra_iterations;

  int keep_iterating = 1;


  // Initialize various parameters and variables:
  errx = 1. ;
  df = f = 1.;
  i_extra = doing_extra = 0;

  if (EXTRA_NEWT_ITER==0) {
    i_extra=-1;
  }
  x_old[0] = x_orig[0] = x[0] ;
  x_old[1] = x_orig[1] = x[1] ;

  vsq = W = W_old = 0.;
  n_iter = 0;

  double xye = con[YE]/con[D];
  double xrho = con[D]/sqrt((1.0/(1.0-x[1])));
  double xtemp = 0.01;
  double xeps = 0.0;
  double xpress = 0.0;
  int keyerr = 0;
  int anyerr = 0;
  int keytemp = 1;

  EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &xrho, &xeps, &xtemp, &xye, &xpress, &keyerr, &anyerr);

  double h_min = xeps + xpress/xrho;

  /* Start the Newton-Raphson iterations : */
  while( keep_iterating ) {
#if DEBUG
    printf("Before func Dx = %g, %g\n", dx[0],dx[1]);
#endif
    (*funcd) (x, dx, resid, jac, &f, &df, lgp, con, temp_guess, tol_x, h_min);  /* returns with new dx, f, df */
#if DEBUG
    printf("After func Dx = %g, %g\n", dx[0],dx[1]);
#endif

    // save old values before calculating the new
    errx = 0.;
    x_old[0] = x[0] ;
    x_old[1] = x[1] ;

    // make the newton step
    x[0] += dx[0]  ;
    x[1] += dx[1]  ;

    // calculate the convergence criterion
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);

    // make sure that the new x[] is physical
    if( x[0] < 0. ) {  x[0] = fabs(x[0]);  }
    else {
     if(x[0] > W_TOO_BIG)  { x[0] = x_old[0] ; }
    }

    if( x[1] < 0. ) {  x[1] = 0.; }
    else {
      if( x[1] >= 1. ) { x[1] = dv; }
    }

    //n_iter++;

    // If we've reached the tolerance level, then just do a few extra iterations
    //  before stopping
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0))
        || (i_extra >= EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }
    n_iter++;

  }   // END of while(keep_iterating)

  // Check for bad untrapped divergences
  if( (!isfinite(f)) ||  (!isfinite(df)) ) {
    c2p_rep->count = n_iter;
    c2p_rep->failed = true;
    return(2);
  }

  if( fabs(errx) <= NEWT_TOL ){
    c2p_rep->count = n_iter;
    c2p_rep->failed = false;
    return(0);
  }
  else if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    c2p_rep->count = n_iter;
    c2p_rep->failed = false;
    return(0);
  }
  else {
    c2p_rep->count = n_iter;
    c2p_rep->failed = true;
    return(1);
  }

}

/*********************************************************************************
   func_vsq():

        -- calculates the residuals, and Newton step for twod_newton_raphson();

     Arguments:
          x   = current value of independent variables (W,vsq) (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/

static void func_vsq(double x[], double dx[], double resid[],
      double jac[][2], double *f, double *df, struct LocGlob *lgp,
      const double * con, double * temp_guess,  const double tol_x, double h_min)
{
  double  W, vsq, p_tmp, dPdvsq, dPdW;
  double t11, t16,t18,t2,t21,   t23,   t24,   t25,   t3,   t35,  t4,   t40;
  double x_new_0, x_new_1, gamma_sq_new, h_new;
  W = x[0];
  vsq = x[1];

  double gamma = 1.0 / sqrt(1.0 - vsq);

  p_tmp = eos_quantities_general2(x, &dPdW, &dPdvsq, lgp, con, temp_guess, tol_x);

  // These expressions were calculated using Mathematica, but made into efficient
  // code using Maple. Since we know the analytic form of the equations, we can
  // explicitly calculate the Newton-Raphson step:

  t2 = -lgp->half_Bsq+dPdvsq;
  t3 = lgp->Bsq+W;
  t4 = t3*t3;
  t23 = 1/W;
  t16 = lgp->QdotBsq*t23*t23;
  t11 = lgp->Qtsq-vsq*t4+t16*(lgp->Bsq+W+W);
  t18 = -lgp->Qdotn-lgp->half_Bsq*(1.0+vsq)+0.5*t16-W+p_tmp;
  t24 = t16*t23;
  t25 = -1.0+dPdW-t24;
  t35 = t25*t3+(lgp->Bsq-2.0*dPdvsq)*(t16+vsq*W)*t23;
  t21 = 1/(t3*t35);
  dx[0] = -(t2*t11+t4*t18)*t21;
  t40 = -2*(vsq+t24)*t3;
  dx[1] = -(-t25*t11+t40*t18)*t21;
  x_new_0 = x[0] + dx[0];
  x_new_1 = x[1] + dx[1];
  gamma_sq_new = 1./(1.-x_new_1);
  double rho_new = con[D]/(sqrt(gamma_sq_new));
  double rho_old = con[D]/gamma;
  h_new = x_new_0/(gamma_sq_new * rho_new) - 1.0;

  jac[0][0] = t40;
  jac[0][1] = -t4;
  jac[1][0] = t25;
  jac[1][1] = t2;
  resid[0] = t11;
  resid[1] = t18;

  *df = -resid[0]*resid[0] - resid[1]*resid[1];

  *f = -0.5 * ( *df );

}


void comp_press_eps_from_x2(double * press, double * eps, const double * x, const double BdotS, const double B_squared, const double * con )
{
  double W = 1./(sqrt(1.0-x[1]));
  double z = x[0];

  *press = -0.5*B_squared/(W*W) - con[TAU] - con[D]+z+B_squared-0.5*pow(BdotS,2)/(z*z);

  *eps   = (z - con[D]*W - (*press)*W*W)/(con[D]*W);

}


void calc_prim_harm(struct LocGlob *lgp, const double * con, double * prim, const double temp_guess, double * x, const double g_cov[4][4], const double g_con[4][4], struct c2p_report * c2p_rep, const double * METRIC, const double * METRIC_PHYS, const double * METRIC_LAP_PSI4, const int iteration, const double xco, const double yco, const double zco){
  DECLARE_CCTK_PARAMETERS;

  /* Recover the primitive variables from the scalars (W,v^2)
   * and conserved variables
   *
   * Eq (23)-(25) in Cerdá-Durán et al. 2008
   */

  double W = 1.0/sqrt(1.0 - x[1]);
  double Z = x[0];

  double B1_cov, B2_cov, B3_cov;
  double con_B1, con_B2, con_B3; 
  con_B1 = con[B1_con]*ONE_OVER_SQRT_4PI; con_B2 = con[B2_con]*ONE_OVER_SQRT_4PI; con_B3 = con[B3_con] * ONE_OVER_SQRT_4PI;
  //con_B1 = con[B1_con]; con_B2 = con[B2_con]; con_B3 = con[B3_con];
  // Lower indices - covariant
  // you can use the 4-metric here because g_{ij} = gamma_{ij} (3-metric)
  B1_cov = g_cov[1][1]*con_B1+g_cov[1][2]*con_B2+g_cov[1][3]*con_B3;
  B2_cov = g_cov[1][2]*con_B1+g_cov[2][2]*con_B2+g_cov[2][3]*con_B3;
  B3_cov = g_cov[1][3]*con_B1+g_cov[2][3]*con_B2+g_cov[3][3]*con_B3;
  //B1_cov = g_cov[1][1]*con[B1_con]+g_cov[1][2]*con[B2_con]+g_cov[1][3]*con[B3_con];
  //B2_cov = g_cov[1][2]*con[B1_con]+g_cov[2][2]*con[B2_con]+g_cov[2][3]*con[B3_con];
  //B3_cov = g_cov[1][3]*con[B1_con]+g_cov[2][3]*con[B2_con]+g_cov[3][3]*con[B3_con];

  prim[RHO] = con[D]/W;
  /*============================ START NEW ============================*/
  //NOTE THE LORENTZ FACTOR! For some reason C2P doesn't include the 
  //Lorentz factor and Illinois does...
  prim[v1_cov] = W*(con[S1_cov] + (lgp->QdotB)*B1_cov/Z)/(Z+lgp->Bsq);
  prim[v2_cov] = W*(con[S2_cov] + (lgp->QdotB)*B2_cov/Z)/(Z+lgp->Bsq);
  prim[v3_cov] = W*(con[S3_cov] + (lgp->QdotB)*B3_cov/Z)/(Z+lgp->Bsq);
  //prim[v1_cov] = (con[S1_cov] + (lgp->QdotB)*B1_cov/Z)/(Z+lgp->Bsq);
  //prim[v2_cov] = (con[S2_cov] + (lgp->QdotB)*B2_cov/Z)/(Z+lgp->Bsq);
  //prim[v3_cov] = (con[S3_cov] + (lgp->QdotB)*B3_cov/Z)/(Z+lgp->Bsq);
  /*============================ END NEW ============================*/

  // calculate P, eps from rho, temp, Ye using EOS
  // in this way P and eps are consistent with the other primitives
  const int keytemp = 0;
  int keyerr = 0;
  int anyerr = 0;

  double xrho, xtemp, xye, xeps, xprs, xentr, xabar;
  xrho = prim[RHO];
  xtemp = temp_guess;
  xye = con[YE]/con[D];
  xeps = 0.0;
  xprs = 0.0;
  xentr = 0.0;
  xabar = 0.0;

  /* using inversion of h to temp (less robust)
  double dpdeps = 0.0;
  double dpdrho = 0.0;
  double xenth = Z/(W*W*xrho); // specific enthalpy h, Z = rho*h*W^2
  EOS_P_from_hrho_dPdrho_dPdeps(xrho, xenth, con, &xtemp, &xeps, &xprs, &dpdrho,
          &dpdeps, &xentr, &xabar, &nEOS_calls);
  if (c2p.evolve_T){
    prim[ENT] = xentr;
    prim[A_BAR] = xabar;
  }
  */

  // not using inversion of h to temp
  comp_press_eps_from_x2(&xprs, &xeps, x, lgp->QdotB, lgp->Bsq, con);
  if (c2p.evolve_T){
    EOS_press_ent_abar(c2p.eoskey,keytemp,xrho,&xeps,&xtemp,xye,&xprs,&xentr,&xabar,&keyerr);
    prim[ENT] = xentr;
    prim[A_BAR] = xabar;
//    prim[MU_HAT] = 0.0;
  } else {
    EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &xrho, &xeps, &xtemp, &xye, &xprs, &keyerr, &anyerr);
  }

  prim[EPS] = xeps;
/*============================ START NEW ============================*/
// different conventions on whether the B-field carries these 4*pi factors
// here we use the Illinois convention
  prim[B1_con] = con[B1_con];
  prim[B2_con] = con[B2_con];
  prim[B3_con] = con[B3_con];
/*============================ END NEW ============================*/
  prim[TEMP] = xtemp;
  prim[PRESS] = xprs;
  prim[WLORENTZ] = W;
  prim[YE] = con[YE]/con[D];

  if( (prim[RHO] <= 0.) || (prim[PRESS] <= 0.) ) {
    c2p_rep->failed=true;
    c2p_rep->c2p_keyerr=5;
  }

    
  if(c2p_rep->c2p_keyerr!=0){c2p_rep->failed=true;}
  else{c2p_rep->failed=false;}

}

/**********************************************************************
  The following routines specify the equation of state. All routines
  above here should be indpendent of the EOS.
**********************************************************************/

/**********************************************************************
  eos_quantities():

      -- returns EOS-related quantities;
 **********************************************************************/

// this is the original routine, only works for ideal gas EOS

// static double eos_quantities(double W, double vsq, double *dpdw, double *dpdvsq, struct LocGlob *lgp)
// {
//   // only works for ideal gas!!!
//   register double ftmp,gtmp;
//
//   double gammaeos = EOS_const[0];
//
//   ftmp = 1. - vsq;
//   gtmp = sqrt(ftmp);
//
//   double gam_m1_o_gam = ((gammaeos-1.)/gammaeos);
//
//
//   *dpdw =  gam_m1_o_gam * ftmp ;    //1)original
//
//   double press = gam_m1_o_gam * ( W * ftmp  -  lgp->Dcon * gtmp );
//   double rho = lgp->Dcon * gtmp;
//   double eps = press / (rho * (gammaeos -1.0));
//
//   *dpdvsq =  gam_m1_o_gam * ( 0.5 * lgp->Dcon/gtmp  -  W ) ;
//
//   return press;  // p
//
// }


static double eos_quantities_general2(const double * x, double *dpdw, double *dpdvsq,
        const struct LocGlob *lgp, const double * con, double * temp_guess,
        const double tol_x)
{

  /*
  - provides p, dp/dW, dp/d(v^2) for given W, v^2, conservatives
  - works for 3-parameter EOS
  - performs an inversion from specific enthalpy to temperature and
    computes the other quantities
  */

  double W   = x[0]; // rho*h*gamma^2
  double vsq = x[1]; // v^2

  double gamma_sq = 1.0/(1. - vsq);
  double gamma = sqrt(gamma_sq);

  double dpdeps = 0.0;
  double dpdrho = 0.0;
  double p = 0.0;
  double eps = 0.0;
  double entr = 0.0;
  double abar = 0.0;
  const double rho = con[D] /gamma;
  double enth = W/gamma_sq/rho; // specific enthalpy h, W = rho*h*gamma^2

#if DEBUG
  printf("***************new iteration****************\n");
  printf("W=%e, gamma=%e, rho=%e, enth=%e\n",W,gamma,rho,enth);
  printf("temp_guess= %g\n", *temp_guess);
#endif

  EOS_P_from_hrho_dPdrho_dPdeps(rho, enth, con, temp_guess, &eps, &p, &dpdrho, &dpdeps, &entr, &abar);

#if DEBUG
  printf("press= %g\n", p);
#endif

  double dpdeps_o_rho= dpdeps/rho;
  *dpdw = ( dpdeps_o_rho/(1.+dpdeps_o_rho) )/gamma_sq;

  double dpdvsq_1 = -0.5*con[D]*gamma*dpdrho;
  double dpdvsq_2 = -0.5*(W + p*gamma_sq)/rho;

  *dpdvsq = (dpdvsq_1 + dpdeps*dpdvsq_2)/(1+dpdeps_o_rho);

  return p;

}
