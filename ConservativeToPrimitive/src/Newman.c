/*@@
   @file      Newman.c
   @date      Sep 17, 2016
   @author    Philipp Moesta, Samantha Chloe Wu
   @desc
   Newman solver
   @enddesc
 @@*/



#include "con2prim.h"


extern struct c2p_steer c2p;


void newman(struct c2p_report * c2p_rep, const double S_squared,
          const double BdotS, const double B_squared, const double * con, double * prim,
	      const double g_cov[4][4], const double g_con[4][4], const double tol_x, double xco, double yco, double zco,
	      int iteration, double * METRIC_LAP_PSI4){


  const int eoskey = c2p.eoskey;
  bool conacc = 0;
  const double sqrtdetg = 1.0;
//  sqrtdetg = sqrt(METRIC_LAP_PSI4[C2P_LAPSE]*METRIC_LAP_PSI4[C2P_PSI6]); // == alpha sqrt{gamma} = alpha Psi^6
  // 2) Compute useful temporary variables
  const double cE = (con[TAU] + con[D])/sqrtdetg;
  double cMsqr = 0.;
  double cT = 0.;
  double cBsqr = 0.;
  double cDetg = sqrtdetg*sqrtdetg;
  double prec = c2p.tol_x;
  cT = BdotS/cDetg; // cT += evolvedVariables.tildeS(i)[s]*evolvedVariables.tildeB(i)[s];
  cMsqr = S_squared/cDetg; //auxiliaryVariables.upperSpatialMetric(i,j)[s]*evolvedVariables.tildeS(i)[s]*evolvedVariables.tildeS(j)[s];
  cBsqr = B_squared;//auxiliaryVariables.lowerSpatialMetric(i,j)[s]*B[i][s]*B[j][s];

  double cP = 0.0;
  double xrho = 0.0;
  double xye = con[YE]/con[D];
  double xeps = 0.0;
  double xtemp = c2p.T_atmo;
//  double xtemp = prim[TEMP];
  int keyerr = 0;
  int anyerr = 0;
  int keytemp = 1;

  xrho = con[D]/prim[WLORENTZ];

  EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &xrho, &xeps, &xtemp, &xye, &cP, &keyerr, &anyerr);

  //Now, begin iterative procedure to derive the primitive variables
  double cPold = cP; // -> setup pressure initial guess
  int step=0;
  double cL;

  const int maxsteps = c2p.max_iterations;
  double AtP[maxsteps] ;    //not length 3 in case of extrap. probs
  double AtR;
  int AtStep=0;
  AtP[0]=cP;

  double cD = 0.5*(cMsqr*cBsqr-cT*cT);
  if(cD<0.) cD=0.;
  bool PfromMaxV = false;

#if DEBUG
 printf("cE = %e, cMsqr = %e, cT = %e, cBsqr = %e, cDetg = %e, cP = %e,con[TAU]=%e,con[D]=%e \n", cE, cMsqr, cT, cBsqr, cDetg, cP, con[TAU], con[D]);
#endif

 do {
    cPold = cP;
    step++;
    //This check is required to ensure that |cos(cPhi)|<1 below
    double cA = cE + cP + cBsqr/2.;
    //Check that the values are physical...
    //What's the right thing to do if cD<0.?
    double cPhi = acos(sqrt(27./4.*cD/cA)/cA);
    double cEps = cA/3.*(1.-2.*cos(2./3.*cPhi+2./3.*M_PI));
    cL = cEps-cBsqr;

    double cVsqr = (cMsqr*cL*cL+cT*cT*(cBsqr+2.*cL))/(cL*(cL+cBsqr)*cL*(cL+cBsqr));
    const double cWsqr = 1./(1.-cVsqr);
    double cH = cL/cWsqr;

    prim[WLORENTZ] = sqrt(cWsqr);
    prim[RHO] = con[D]/prim[WLORENTZ]; //    rho[s] = tildeD[s]/(sqrtDetg[s]*W[s]);

    // We'll need a 1D root solve to find the temperature from density
    // and enthalpy.
    const int keytemp = 0;
    const double precEOS = c2p.eos_prec;
    double xeps = 0.0;
    double xye = prim[YE];
    double xprs = 0.0;
    double xrho = prim[RHO];
    int anyerr = 0;
    int keyerr = 0;
    const int npoints = 1;

 // define other dummy variables as needed for input arguments

#if DEBUG
    printf("precall cH = %g\n",cH);
#endif


    EOS_Omni_press_from_rhoenthalpy(eoskey, keytemp, precEOS, 1, &xrho, &xeps, &xtemp, &xye, &xprs, &cH, &keyerr, &anyerr);
//    EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &xrho, &xeps, &xtemp, &xye, &xprs, &keyerr, &anyerr);

    cP = xprs;
    prim[EPS] = xeps;
    prim[TEMP] = xtemp;

#if DEBUG
    printf("POS xprs= %e,cPhi = %e, cEps = %e, cH = %e, cVsqr = %e, cL = %e, cP = %e,prim[WLORENTZ]=%e, prim[RHO]=%e xeps=%e xtemp=%e\n", xprs, cPhi, cEps, cH, cVsqr, cL, cP, prim[WLORENTZ], prim[RHO],xeps,xtemp);
#endif

#if DEBUG
    printf("At step %i, P = %0.12e, dP = %e, cP+cPold = %e, threshold = %e \n",step,cP,fabs(cP-cPold),fabs(cP+cPold),fabs(cP-cPold)/(cP+cPold));
#endif

    AtStep++;
    AtP[AtStep]=cP;

    if(AtStep>=2) {   //Aitken extrapolation

       AtR = (AtP[AtStep]-AtP[AtStep-1])/(AtP[AtStep-1]-AtP[AtStep-2]);
#if DEBUG
    printf("At step %i, AtP[0] = %e, AtP[1] = %e, AtP[2]=%e, AtR=%e \n",step,AtP[0], AtP[1], AtP[2],AtR);
#endif
	if(AtR<1. && AtR>0.) {
        cP=AtP[AtStep-1]+(AtP[AtStep]-AtP[AtStep-1])/(1.-AtR);
#if DEBUG
        printf("At step in acceleration  %i, P = %e, dP = %e \n",AtStep,cP,fabs(cP-cPold));
#endif
        AtStep=0;
        conacc = 1;
        AtP[0]=cP;   //starting value for next Aitken extrapolation
      }
    }
  }
  while(fabs(cP-cPold)>prec*(cP+cPold) && step<maxsteps);

  if (conacc==1) {     //converged on an extrap. so recompute vars
    const double cA = cE + cP + cBsqr/2.;
    const double cPhi = acos(sqrt(27./4.*cD/cA)/cA);
    const double cEps = cA/3.*(1.-2.*cos(2./3.*cPhi+2./3.*M_PI));
    cL = cEps-cBsqr;
    const double cVsqr =
      (cMsqr*cL*cL+cT*cT*(cBsqr+2.*cL))/(cL*(cL+cBsqr)*cL*(cL+cBsqr));
    if(cVsqr<0. || cVsqr>=1.) {
    }

    const double cWsqr = 1./(1.-cVsqr);

    prim[WLORENTZ] = sqrt(cWsqr); //    W[s] = sqrt(cWsqr);
    prim[RHO] = con[D]/prim[WLORENTZ]; //    rho[s] = tildeD[s]/(sqrtDetg[s]*W[s]);
    prim[PRESS]=cP;

    prim[B1_con] = con[B1_con];
    prim[B2_con] = con[B2_con];
    prim[B3_con] = con[B3_con];

    const int keytemp = 0;
    const double prec = 0.0;
    double xeps = 0.0;
    double xtemp = 0.0;
    double xye = 0.0;
    double xprs = 0.0;
    int  anyerr = 0;
    int keyerr = 0;
    double h = 0.0;
  }

 //Compute v^i
  const double cS = cT/cL;

  prim[B1_con] = con[B1_con];
  prim[B2_con] = con[B2_con];
  prim[B3_con] = con[B3_con];

  // Lower indices - covariant

  const double B1_cov = g_cov[1][1]*con[B1_con]+g_cov[1][2]*con[B2_con]+g_cov[1][3]*con[B3_con];
  const double B2_cov = g_cov[1][2]*con[B1_con]+g_cov[2][2]*con[B2_con]+g_cov[2][3]*con[B3_con];
  const double B3_cov = g_cov[1][3]*con[B1_con]+g_cov[2][3]*con[B2_con]+g_cov[3][3]*con[B3_con];

  prim[v1_cov] = prim[WLORENTZ]*(cS * B1_cov + con[S1_cov]) / (sqrtdetg*(cL+cBsqr)); // v[i][s] = cS*evolvedVariables.tildeB(i)[s];
  prim[v2_cov] = prim[WLORENTZ]*(cS * B2_cov + con[S2_cov]) / (sqrtdetg*(cL+cBsqr)); // v[i][s] = cS*evolvedVariables.tildeB(i)[s];
  prim[v3_cov] = prim[WLORENTZ]*(cS * B3_cov + con[S3_cov]) / (sqrtdetg*(cL+cBsqr)); // v[i][s] = cS*evolvedVariables.tildeB(i)[s];

  prim[PRESS]=cP;
  prim[YE] = con[YE]/con[D];

  if(c2p.eoskey==4){
     CCTK_REAL xrhoS=prim[RHO]; CCTK_REAL xepsS=prim[EPS];     CCTK_REAL; CCTK_REAL xtempS = prim[TEMP];
     CCTK_REAL xyeS=prim[YE];   CCTK_REAL xpressS=prim[PRESS]; CCTK_REAL; CCTK_REAL xentS = 0.0; CCTK_INT keytempS = 1;
     CCTK_REAL xcs2S=0.0; CCTK_REAL xdedtS = 0.0; CCTK_REAL xdpderhoS = 0.0; CCTK_REAL xdpdrhoeS = 0.0; CCTK_REAL xmunuS = 0.0;
     EOS_Omni_short(c2p.eoskey, keytempS, c2p.eos_prec, 1, &xrhoS, &xepsS, &xtempS, &xyeS, &xpressS, &xentS, &xcs2S, &xdedtS, &xdpderhoS, &xdpdrhoeS, &xmunuS, &keyerr, &anyerr);
  prim[ENT]      = xentS;
  prim[MU_HAT]   = 0.0;
  }

  c2p_rep->count = step;
  if (step >= maxsteps) c2p_rep->failed = true;

}
