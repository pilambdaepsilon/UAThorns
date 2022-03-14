#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "con2prim.h"


CCTK_REAL P_from_Prim(CCTK_REAL const * prim){

    // Calculate pressure from primitive variables, given rho, T, Ye

    CCTK_INT keyerr = 0;
    CCTK_INT anyerr = 0;
    const CCTK_INT keytemp = 1;

    CCTK_REAL xeps,xprs,xrho,xtemp,xye;
    xrho = prim[RHO];
    xtemp = prim[TEMP];
    xeps = 0.0;
    xye = prim[YE];
    xprs = 0.0;

#if DEBUG
    printf("P_prim: Before call: Temp: %g Rho: %g Eps: %g Ye: %g p: %g\n",xtemp,xrho,xeps,xye,xprs);
#endif

    EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &xrho, &xeps, &xtemp, &xye, &xprs, &keyerr, &anyerr);
    printf("Pressure in P from prim %g \n", xprs);

#if DEBUG
    if (keyerr != 0) printf("P_prim keyerr: %d, Temp: %g Rho: %g Eps: %g Ye: %g p: %g\n",keyerr,xtemp,xrho,xeps,xye,xprs);
#endif

    printf("This is the function that has been called, in NR_utils\n");

    return xprs;
}

CCTK_REAL T_from_Prim(CCTK_REAL const * prim){

    // Calculate temperature from primitive variables, given rho, eps, Ye
    // only required by standalone version

    CCTK_INT keyerr = 0;
    CCTK_INT anyerr = 0;
    const CCTK_INT keytemp = 0;

    CCTK_REAL xeps,xprs,xrho,xtemp,xye;
    xrho = prim[RHO];
    xtemp = 0.0;
    xeps = prim[EPS];
    xye = prim[YE];
    xprs = 0.0;

    EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &xrho, &xeps, &xtemp, &xye, &xprs, &keyerr, &anyerr);

#if DEBUG
    if (keyerr != 0) printf("T_prim keyerr: %d, Temp: %g Rho: %g Eps: %g Ye: %g p: %g\n",keyerr,xtemp,xrho,xeps,xye,xprs);
#endif

    return xtemp;
}

CCTK_REAL eps_from_Prim(CCTK_REAL const * prim){

    // Calculate specific internal energy from primitive variables, given rho, eps, Ye
    // only required by standalone version

    CCTK_INT keyerr = 0;
    CCTK_INT anyerr = 0;
    const CCTK_INT keytemp = 1;

    CCTK_REAL xeps,xprs,xrho,xtemp,xye;
    xrho = prim[RHO];
    xtemp = prim[TEMP];
    xeps = 0.0;
    xye = prim[YE];
    xprs = 0.0;

#if DEBUG
    printf("eps_prim: Before call: Temp: %g Rho: %g Eps: %g Ye: %g p: %g\n",xtemp,xrho,xeps,xye,xprs);
#endif

    EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &xrho, &xeps, &xtemp, &xye, &xprs, &keyerr, &anyerr);

    printf("esp in eps from prim %g \n", xeps);

#if DEBUG
    if (keyerr != 0) printf("eps_prim keyerr: %d, Temp: %g Rho: %g Eps: %g Ye: %g p: %g\n",keyerr,xtemp,xrho,xeps,xye,xprs);
#endif

    return xeps;
}

void EOS_press_ent_abar(int eoskey, int keytemp, double rho, double * eps,
        double * temp, double ye, double * prs, double * ent, double * abar,
        int * keyerr){

  const CCTK_INT eos_key = eoskey;
  const CCTK_INT havetemp = keytemp;
  const CCTK_REAL prec = c2p.eos_prec;
  const CCTK_INT npoints = 1;
  CCTK_INT keyerr_loc;
  CCTK_INT anyerr = 0;
  CCTK_REAL xrho, xeps, xtemp, xye, xprs, xent, xabar, xcs2, xdedt, xdpderho, xdpdrhoe,
    xxa, xxh, xxn, xxp, xzbar, xmue, xmun, xmup, xmuhat;

  keyerr_loc = 0;

  xrho = rho;
  xtemp = *temp;
  xeps = *eps;
  xye = ye;
  xprs = 0.0;
  xent = 0.0;
  xabar = 0.0;

  EOS_Omni_full(eos_key,havetemp,prec,npoints,&xrho,&xeps,&xtemp,&xye,&xprs,&xent,
                &xcs2, &xdedt, &xdpderho, &xdpdrhoe, &xxa, &xxh, &xxn, &xxp, &xabar,
                &xzbar, &xmue, &xmun, &xmup, &xmuhat, &keyerr_loc, &anyerr);

  *eps  = xeps;
  *temp = xtemp;
  *prs  = xprs;
  *ent  = xent;
  *abar = xabar;
  *keyerr = keyerr_loc;
}


void EOS_Press_Cold(const double* rho, const double *ye, double * prs, double * xeps, int * keyerr, int *anyerr){

  int npoints  = 1;
  int havetemp = 1;
  double temp  = c2p.eos_temp_min;

  EOS_Omni_press(c2p.eoskey, havetemp, c2p.eos_prec, npoints, rho, xeps, &temp, ye, prs, keyerr, anyerr);

}


void EOS_EP_dEdr_dEdt_dPdr_dPdt(double * x, const double * con,
          double * Eprim, double * Pprim, double * dEdrho, double * dEdt, double * dPdrho,
          double * dPdt, int stepsize){

  /*  Compute partial derivatives of specific internal energy and pressure with respect
   *  to density and temperature, based on primitives computed from Newton-Raphson state
   *  vector x and conservatives
   */

  const double W = x[0];
  const double Z = x[1];
  const double T = x[2];

  const CCTK_INT eos_key = c2p.eoskey;
  const CCTK_INT keytemp = 1;
  const CCTK_INT npoints = 1;
  const CCTK_INT step_size = stepsize;
  const CCTK_REAL prec = c2p.eos_prec;
  CCTK_INT keyerr;
  CCTK_INT anyerr;
  CCTK_REAL xrho,xtemp,xeps,xye,xprs,xdedrho,xdpdrho,xdedt,xdpdt;

  keyerr = 0;

  xrho = con[D]/W;
  xtemp = T;
  xeps = 0.0;
  xye = con[YE]/con[D];
  xprs = 0.0;

  xdedrho = 0.0;
  xdpdrho = 0.0;
  xdedt = 0.0;
  xdpdt = 0.0;


  EOS_Omni_dpdrho_dpdt_dedrho_dedt(eos_key, keytemp, c2p.eos_prec, npoints, \
                                   &xrho, &xtemp, &xye, &xeps, &xprs, &xdpdrho, &xdpdt, &xdedrho, &xdedt, &keyerr, &anyerr );

  printf("EOS_EP_dEdr_dEdt_dPdr_dPdt: Temp: %g Rho: %g Eps: %g Ye: %g p: %g\n",xtemp,xrho,xeps,xye,xprs);

#if DEBUG
  if (keyerr != 0) printf("EOS_EP_dEdr_dEdt_dPdr_dPdt: keyerr: %d, Temp: %g Rho: %g Eps: %g Ye: %g p: %g\n",keyerr,xtemp,xrho,xeps,xye,xprs);
#endif

  *dEdrho = xdedrho;
  *dEdt   = xdedt;
  *dPdrho = xdpdrho;
  *dPdt   = xdpdt;
  *Eprim  = xeps;
  *Pprim  = xprs;

}


void EOS_P_from_hrho_dPdrho_dPdeps(const double rho, const double enth, const double * con,
          double * temp_guess, double * eps, double * press, double * dPdrho,
          double * dPdeps, double * entr, double * abar){

  /*  Perform inversion from specific enthalpy to temperature and compute partial derivatives
   *  of pressure wrt density and specific internal energy
   */

  const CCTK_INT eos_key = c2p.eoskey;
  const CCTK_INT step_size = 1;
  const CCTK_REAL prec = c2p.eos_prec;
  const CCTK_REAL xrho = rho;
  const CCTK_REAL xenth = enth;
  const CCTK_REAL xye = con[YE]/con[D];
  CCTK_INT keyerr;
  CCTK_INT anyerr;
  CCTK_INT npoints = 1;
  CCTK_REAL xtemp,xeps,xprs,xdpdrho,xdpdeps,xentr,xabar;

  keyerr = 0;
  xtemp = *temp_guess;
  xeps = 0.0;
  xprs = 0.0;
  xdpdeps = 0.0;
  xdpdrho = 0.0;
  xentr = 0.0;
  xabar = 0.0;

  EOS_Omni_dpdrho_dpdeps_hinv(eos_key, prec, npoints, &xrho, &xenth, &xtemp, &xye, &xeps, &xprs, &xdpdrho,
                              &xdpdeps, &keyerr, &anyerr);

#if DEBUG
  if (keyerr != 0) printf("EOS_P_from_hrho_dPdr_dPdeps: keyerr: %d, Temp_guess: %g Rho: %g enthalpy: %g Ye: %g p: %g\n",keyerr,xtemp,xrho,xenth,xye,xprs);
#endif

  *temp_guess = xtemp;
  *dPdeps = xdpdeps;
  *dPdrho = xdpdrho;
  *press  = xprs;
  *eps = xeps;
  *entr = xentr;
  *abar = xabar;

}


void EOS_EP_dEdr_dEdt_dPdr_dPdt_2D(const double rho2D, const double temp2D, const double * con,
          double * Eprim, double * Pprim, double * dEdrho, double * dEdt, double * dPdrho,
          double * dPdt, int stepsize){

  /*  Compute partial derivatives of specific internal energy and pressure with respect
   *  to density and temperature, based on primitives computed from Newton-Raphson state
   *  vector x and conservatives
   */

  const CCTK_INT eos_key = c2p.eoskey;
  const CCTK_INT keytemp = 1;
  const CCTK_INT step_size = stepsize;
  const CCTK_REAL prec = c2p.eos_prec;
  CCTK_INT keyerr;
  CCTK_INT anyerr;
  const CCTK_INT npoints = 1;
  CCTK_REAL xrho,xtemp,xeps,xye,xprs,xdedrho,xdpdrho,xdedt,xdpdt;

  keyerr = 0;

  xrho = rho2D; //con[D]/x[0];
  xtemp = temp2D; //x[1];
  xeps = 0.0;
  xye = con[YE]/con[D];
  xprs = 0.0;

  xdedrho = 0.0;
  xdpdrho = 0.0;
  xdedt = 0.0;
  xdpdt = 0.0;

  EOS_Omni_dpdrho_dpdt_dedrho_dedt(eos_key, keytemp, c2p.eos_prec, npoints, \
                                   &xrho, &xtemp, &xye, &xeps, &xprs, &xdpdrho, &xdpdt, &xdedrho, &xdedt, &keyerr, &anyerr );

#if DEBUG
  if (keyerr != 0) printf("EOS_EP_dEdr_dEdt_dPdr_dPdt: keyerr: %d, Temp: %g Rho: %g Eps: %g Ye: %g p: %g\n",keyerr,xtemp,xrho,xeps,xye,xprs);
#endif

  *dEdrho = xdedrho;
  *dEdt   = xdedt;
  *dPdrho = xdpdrho;
  *dPdt   = xdpdt;
  *Eprim  = xeps;
  *Pprim  = xprs;

}
