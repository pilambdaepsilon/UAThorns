/*@@
   @file      NR_utils.c
   @date      Sep 17, 2016
   @author    Daniel Siegel, Philipp Moesta
   @desc
   Utils for NR.c
   @enddesc
 @@*/

#include <assert.h>
#include "con2prim.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

extern struct c2p_steer c2p;

void EOSprim_from_rhoTYe(const int keytemp, double * prim){

    // Calculate EOS primitives, given rho, T, Ye

    assert(keytemp < 2);
    int keyerr = 0;
    int anyerr = 0;

    double xeps,xprs,xrho,xtemp,xye;
    xrho = prim[RHO];
    xye = prim[YE];
    xprs = 0.0;

    if (c2p.evolve_T && c2p.evolve_Ye) {
      double xabar = 0.0;
      double xent = 0.0;
      if (keytemp == 1) {
        xtemp = prim[TEMP];
        xeps = 0.0;
      } else {
        xtemp = 0.0;
        xeps = prim[EPS];
      }

      EOS_press_ent_abar(c2p.eoskey,keytemp,xrho,&xeps,&xtemp,xye, &xprs,&xent,&xabar,&keyerr);

      prim[TEMP] = xtemp;
      prim[EPS] = xeps;
      prim[ENT] = xent;
      prim[A_BAR] = xabar;
      prim[MU_HAT] = 0.0;

    } else {
      xtemp = 0.0;
      xeps = prim[EPS];
      assert(keytemp == 0);
      EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &xrho, &xeps, &xtemp, &xye, &xprs, &keyerr, &anyerr);


      prim[TEMP] = xtemp;
    }

    prim[PRESS] = xprs;

    if (keyerr != 0) printf("EOSprim_from_rhoTYe keyerr: %d, Temp: %g Rho: %g Eps: %g Ye: %g p: %g\n",keyerr,xtemp,xrho,xeps,xye,xprs);

}


void EP_dEdW_dEdZ_dEdT(double * Eprim, double * Pprim, double * dEdW, double * dEdZ, double * dEdT,
       double * dpdrho, double * dpdT, double * x, const double * con, int stepsize){

  // Compute partial derivatives of specific internal energy and pressure with respect
  // to W and Z
  // Note: E = eps - eps_EOS

  double W = x[0];
  double Z = x[1];
  double T = x[2];

  double epsEOS,pEOS,depsEOSdrho,dpEOSdrho,depsEOSdt,dpEOSdt;
  epsEOS = 0.0;
  pEOS = 0.0;
  depsEOSdrho = 0.0;
  dpEOSdrho=0.0;
  dpEOSdt = 0.0;

  /* need partial derivatives of specific internal energy and pressure wrt density and
   * temperature. Those need to be based on primitives computed from Newton-Raphson state
   * vector x and conservatives
   */
  EOS_EP_dEdr_dEdt_dPdr_dPdt(x,con,&epsEOS,&pEOS,&depsEOSdrho,&depsEOSdt,&dpEOSdrho,&dpEOSdt,stepsize);

  // Further partial derivatives
  double drhodW = -con[D]/(W*W);
  double depsdW = -Z/(con[D]*W*W)-pEOS/con[D] + dpEOSdrho/W;
  double depsdP = -W/con[D];
  double drhodZ = 0.0;
  double depsdZ = 1.0/(con[D]*W);

  *Eprim = epsEOS;
  *Pprim = pEOS;
  *dEdW = depsdW - depsEOSdrho*drhodW;
  *dEdZ = depsdZ;
  *dEdT = depsdP*dpEOSdt - depsEOSdt;
  *dpdrho = dpEOSdrho;
  *dpdT = dpEOSdt;
}


void EP_dPdW_dPdZ_dPdT(double * Eprim, double * Pprim, double * dPdW, double * dPdZ, double * dPdT,
       double * dpdrho, double * dpdT, double * x, const double * con, int stepsize){

  // Compute partial derivatives of specific internal energy and pressure with respect
  // to W and Z

  double W = x[0];
  double Z = x[1];
  double T = x[2];

  double epsEOS,pEOS,depsEOSdrho,dpEOSdrho,depsEOSdt,dpEOSdt;
  epsEOS = 0.0;
  pEOS = 0.0;
  depsEOSdrho = 0.0;
  dpEOSdrho=0.0;
  dpEOSdt = 0.0;

  /* need partial derivatives of specific internal energy and pressure wrt density and
   * temperature. Those need to be based on primitives computed from Newton-Raphson state
   * vector x and conservatives
   */
  EOS_EP_dEdr_dEdt_dPdr_dPdt(x,con,&epsEOS,&pEOS,&depsEOSdrho,&depsEOSdt,&dpEOSdrho,&dpEOSdt,stepsize);

  // Further partial derivatives
  double drhodW = -con[D]/(W*W);
  double drhodZ = 0.0;
  double dpdW = (con[D]*(1.0 + epsEOS + con[D]*depsEOSdrho/W) - 2.0*Z/W) / (W*W);
  double dpdeps = -con[D]/W;
  double dpdZ = 1.0/(W*W);

  *Eprim = epsEOS;
  *Pprim = pEOS;
  *dPdW =  dpdW - dpEOSdrho*drhodW;
  *dPdZ = dpdZ;
  *dPdT = dpdeps*depsEOSdt - dpEOSdt;
  *dpdrho = dpEOSdrho;
  *dpdT = dpEOSdt;
}
