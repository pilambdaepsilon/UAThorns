#include "con2prim.h"
#include "brent.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

extern struct c2p_steer c2p;
extern struct brent_report breport;


void calc_prim(const double x, const double * con, const double * param, const double temp_guess, double * prim,
      const double S_squared, const double BdotS,
      const double B_squared, const double g_cov[4][4], struct c2p_report * c2p_rep,const int index,double xco, double yco, double zco,int iteration){

  // Recover the primitive variables prim from x, q, r, s, t, con

  double q = param[par_q];
  double r = param[par_r];
  double s = param[par_s];
  double t = param[par_t];

  double Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t ) / ( x*x*(x+s)*(x+s) );
  Wminus2 = fmin(fmax(Wminus2 ,1e-10 ), 1.0-1e-10);
  double W= pow(Wminus2, -0.5);

  double rho = con[D]/W;

  double eps = W - 1.0 + (1.0-W*W)*x/W + W*(q - s + t*t/(2*x*x) + s/(2*W*W)  );
  eps=fmax(eps, c2p.eos_eps_min);

  double ye = con[YE]/con[D];

  int keytemp = 0;
  int keyerr = 0;
  int anyerr = 0;
  int nEOScalls=0;
  double temp = temp_guess;
  double press = 0.0;

  if (c2p.evolve_T){
    double ent = 0.0;
    double abar = 0.0;
    EOS_press_ent_abar(c2p.eoskey,keytemp,rho,&eps,&temp,ye,&press,&ent,&abar,&keyerr);

    c2p_rep->nEOScalls += nEOScalls;
    prim[ENT] = ent;
    prim[A_BAR] = abar;
    prim[MU_HAT] = 0.0;
  } else {
    EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &rho, &eps, &temp, &ye, &press, &keyerr, &anyerr);
    c2p_rep->nEOScalls += nEOScalls;
  }

  double B1_cov, B2_cov, B3_cov;
  const double Z= x*rho*W;
  double con_B1, con_B2, con_B3; 
  con_B1 = con[B1_con]*ONE_OVER_SQRT_4PI; con_B2 = con[B2_con]*ONE_OVER_SQRT_4PI; con_B3 = con[B3_con] * ONE_OVER_SQRT_4PI;
  // Lower indices - covariant
  B1_cov = g_cov[1][1]*con_B1+g_cov[1][2]*con_B2+g_cov[1][3]*con_B3;
  B2_cov = g_cov[1][2]*con_B1+g_cov[2][2]*con_B2+g_cov[2][3]*con_B3;
  B3_cov = g_cov[1][3]*con_B1+g_cov[2][3]*con_B2+g_cov[3][3]*con_B3;
  //B1_cov = g_cov[1][1]*con[B1_con]+g_cov[1][2]*con[B2_con]+g_cov[1][3]*con[B3_con];
  //B2_cov = g_cov[1][2]*con[B1_con]+g_cov[2][2]*con[B2_con]+g_cov[2][3]*con[B3_con];
  //B3_cov = g_cov[1][3]*con[B1_con]+g_cov[2][3]*con[B2_con]+g_cov[3][3]*con[B3_con];

  prim[RHO] = rho;
  prim[v1_cov] = W*(con[S1_cov] + (BdotS)*B1_cov/Z)/(Z+B_squared);
  prim[v2_cov] = W*(con[S2_cov] + (BdotS)*B2_cov/Z)/(Z+B_squared);
  prim[v3_cov] = W*(con[S3_cov] + (BdotS)*B3_cov/Z)/(Z+B_squared);
  prim[EPS] = eps;
  prim[B1_con] = con[B1_con];
  prim[B2_con] = con[B2_con];
  prim[B3_con] = con[B3_con];
  prim[TEMP] = temp;
  prim[YE] = ye;
  prim[PRESS] = press;
  prim[WLORENTZ] = W;

  if(c2p.eoskey==4){
     CCTK_REAL xrhoS=prim[RHO]; CCTK_REAL xepsS=prim[EPS]; CCTK_REAL xtempS = prim[TEMP]; CCTK_REAL xyeS=prim[YE]; 
     CCTK_REAL xpressS=prim[PRESS]; CCTK_REAL; CCTK_REAL xentS = 0.0; CCTK_REAL xabarS = 0.0; 
     CCTK_REAL xcs2S=0.0; CCTK_REAL xdedtS = 0.0; CCTK_REAL xdpderhoS = 0.0; CCTK_REAL xdpdrhoeS = 0.0; CCTK_REAL xmunuS = 0.0;
     CCTK_REAL xxaS=0.0; CCTK_REAL xxhS = 0.0; CCTK_REAL xxnS = 0.0; CCTK_REAL xxpS = 0.0; CCTK_REAL xzbarS = 0.0;
     CCTK_REAL xmueS=0.0; CCTK_REAL xmunS = 0.0; CCTK_REAL xmupS = 0.0; CCTK_REAL xmuhatS = 0.0;

     CCTK_INT keytempS = 1;

     //TODO: Find a stable way of getting muhat. EOS_Omni_full sometimes interpolates to NaNs in the chemical potential muhat
     EOS_Omni_short(c2p.eoskey, keytempS, c2p.eos_prec, 1, &xrhoS, &xepsS, &xtempS, &xyeS, &xpressS, &xentS, &xcs2S, &xdedtS, &xdpderhoS, &xdpdrhoeS, &xmunuS, &keyerr, &anyerr);
//     EOS_Omni_full(c2p.eoskey, keytempS, c2p.eos_prec, 1, &xrhoS, &xepsS, &xtempS, &xyeS, &xpressS, &xentS, &xcs2S, &xdedtS, &xdpderhoS, &xdpdrhoeS, &xxaS, &xxhS, &xxnS, &xxpS, 
//		      &xabarS, &xzbarS, &xmueS, &xmunS, &xmupS, &xmuhatS, &keyerr, &anyerr);
     prim[ENT]      = xentS;
     prim[MU_HAT]   = 0.0;
//     prim[MU_HAT]   = fmin(xmuhatS, 300);
  }
  else{
     prim[ENT]      = 1.0;
     prim[MU_HAT]   = 0.0;
  }

  double cBsqr=B_squared; double cS = BdotS; double cVsqr = 1.- 1./(W*W);

}


double func_root(double x, double * param, struct c2p_report * c2p_rep, bool use_epsmin, double * temp_guess){

  // computes f(x) from x and q,r,s,t

  int keytemp = 0;
  int keyerr = 0;
  int anyerr = 0;

  double P = 0.0;
  double temp = *temp_guess;

  const double ye = param[conYE]/param[conD];
  const double q= param[par_q];
  const double r= param[par_r];
  const double s= param[par_s];
  const double t= param[par_t];

  // (i)
  double Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t  )/ (x*x*(x+s)*(x+s));

  Wminus2 = fmin(fmax(Wminus2 ,1e-10 ), 1-1e-10);
  const double W= pow(Wminus2, -0.5);

  // (ii)
  double rho = param[conD]/W;

  // (iii)
  double eps;
  eps = W - 1.0 + (1.0-W*W)*x/W + W*(q - s + t*t/(2*x*x) + s/(2*W*W));
  eps = fmax(eps, c2p.eos_eps_min);

  // (iv)
  int nEOScalls=0;
  EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &rho, &eps, &temp, &ye, &P, &keyerr, &anyerr);
  c2p_rep->nEOScalls += nEOScalls;


  double ans = x- (1.0 + eps + P/rho)*W;

  return ans;
}



void palenzuela(struct c2p_report * c2p_rep, const double S_squared, const double BdotS,
      const double B_squared, const double * con, double * prim, const double g_con[4][4],
      const double g_cov[4][4], const double tol_x, bool use_epsmin, const int index, double xco, double yco, double zco, int iteration)
{

  // main root finding routine
  // using scheme of Palenzuela et al. 2015, PRD, 92, 044045

  // some quantities computed from the conservatives,
  // defined as order unity quantities
  double q= con[TAU]/con[D];
  double r= S_squared/(con[D]*con[D]);
  double s= B_squared/con[D];
  double t= BdotS/(pow(con[D],1.5));

  double param[6];
  param[par_q] = q;
  param[par_r] = r;
  param[par_s] = s;
  param[par_t] = t;
  param[conD]  = con[D];
  param[conYE] = con[YE];

  // bracket for x
  double xlow = 1.0+q-s;
  double xup = 2.0+2.0*q-s;

  // initial guess for temperature
  double temp_guess=prim[TEMP];

  // reset #EOS calls counter
  c2p_rep->nEOScalls = 0;

  // find x, this is the recovery process
  double x = zbrent(*func_root, param, &temp_guess, xlow, xup, tol_x, c2p_rep, use_epsmin);

  // calculate final set of primitives
  calc_prim(x, con, param, temp_guess, prim, S_squared, BdotS, B_squared, g_cov, c2p_rep,index,xco,yco,zco,iteration);


}
