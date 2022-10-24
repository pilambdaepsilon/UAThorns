/*@@
   @file      con2prim_MHD.c
   @date      Sep 17, 2016
   @author    Daniel Siegel, Philipp Moesta
   @desc
   This file provides the main con2prim routine con2prim_MHD_ along with some utilities,
   which handle the actual recovery process and error policy. Note that this code is
   independent of the CACTUS infrastructure and does neither assume a specific con2prim
   scheme nor a specific EOS. Presently, the code calls a recovery routine based on
   3D Newton-Raphson root finding, which implements a scheme based on
   Cerda-Duran et al. (2008) that should work with any EOS given in terms of rho, T, Ye.
   However, any other recovery scheme can be used, in principle. Note that this code
   can be directly tested with the attached standalone wrapper by setting STANDALONE 1
   in the header file con2prim.h.
   @enddesc
 @@*/


#include "con2prim.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#define MAXNUMVARS 20
#define SQR(x) ((x) * (x))
extern struct c2p_steer c2p;

void calc_B2_S2_BdotS(const double * con, const double g_con[4][4],
      const double g_cov[4][4], double * B_squared, double * S_squared, double * St_squared, double * BdotS, double * METRIC, double * METRIC_PHYS, double * METRIC_LAP_PSI4)
{

  double B1_cov, B2_cov, B3_cov, S1_con, S2_con, S3_con;
  double con_B1, con_B2, con_B3; 

  con_B1 = con[B1_con]*ONE_OVER_SQRT_4PI; con_B2 = con[B2_con]*ONE_OVER_SQRT_4PI; con_B3 = con[B3_con] * ONE_OVER_SQRT_4PI;
  // Lower indices - covariant
  B1_cov = g_cov[1][1]*con_B1+g_cov[1][2]*con_B2+g_cov[1][3]*con_B3;
  B2_cov = g_cov[1][2]*con_B1+g_cov[2][2]*con_B2+g_cov[2][3]*con_B3;
  B3_cov = g_cov[1][3]*con_B1+g_cov[2][3]*con_B2+g_cov[3][3]*con_B3;
  //B1_cov = g_cov[1][1]*con[B1_con]+g_cov[1][2]*con[B2_con]+g_cov[1][3]*con[B3_con];
  //B2_cov = g_cov[1][2]*con[B1_con]+g_cov[2][2]*con[B2_con]+g_cov[2][3]*con[B3_con];
  //B3_cov = g_cov[1][3]*con[B1_con]+g_cov[2][3]*con[B2_con]+g_cov[3][3]*con[B3_con];

  // Raise indices - contravariant -> we need to use the 3-metric for this...
  // it's ok for the covariant 4-metric but g^{ij} = gamma^{ij} - beta^{i}beta^{j}/alpha^2
  // where gamma^{ij} is the 3-metric
  S1_con = METRIC_PHYS[C2P_GUPXX]*con[S1_cov]+METRIC_PHYS[C2P_GUPXY]*con[S2_cov]+METRIC_PHYS[C2P_GUPXZ]*con[S3_cov];
  S2_con = METRIC_PHYS[C2P_GUPXY]*con[S1_cov]+METRIC_PHYS[C2P_GUPYY]*con[S2_cov]+METRIC_PHYS[C2P_GUPYZ]*con[S3_cov];
  S3_con = METRIC_PHYS[C2P_GUPXZ]*con[S1_cov]+METRIC_PHYS[C2P_GUPYZ]*con[S2_cov]+METRIC_PHYS[C2P_GUPZZ]*con[S3_cov];
  //S1_con = g_con[1][1]*con[S1_cov]+g_con[1][2]*con[S2_cov]+g_con[1][3]*con[S3_cov];
  //S2_con = g_con[1][2]*con[S1_cov]+g_con[2][2]*con[S2_cov]+g_con[2][3]*con[S3_cov];
  //S3_con = g_con[1][3]*con[S1_cov]+g_con[2][3]*con[S2_cov]+g_con[3][3]*con[S3_cov];

  // Need to calculate for (21) and (22) in Cerda-Duran 2008
  // B * S = B^i * S_i -> this is the same as QdotB, up to that extra factor of 1/sqrt(4*pi)
  *BdotS = (con_B1*con[S1_cov] + con_B2*con[S2_cov] + con_B3*con[S3_cov]);
  //*BdotS = con[B1_con]*con[S1_cov] + con[B2_con]*con[S2_cov] + con[B3_con]*con[S3_cov];

  // B^2 = B^i * B_i -> for consistency with Illinois, include this 1/(4*pi) factor
  *B_squared = (con_B1*B1_cov + con_B2*B2_cov + con_B3*B3_cov);
  //*B_squared = con[B1_con]*B1_cov + con[B2_con]*B2_cov + con[B3_con]*B3_cov;

  // S^2 = S^i * S_i -> 
  *S_squared = S1_con*con[S1_cov] + S2_con*con[S2_cov] + S3_con*con[S3_cov];
  //*S_squared = S1_con*con[S1_cov] + S2_con*con[S2_cov] + S3_con*con[S3_cov];

  *St_squared = *S_squared;
}


void reset_to_atmosphere(double * prim, const double * con)
{

  // reset primitives to atmosphere values

  prim[RHO] = c2p.rho_atmo;
  prim[v1_cov] = 0.0;
  prim[v2_cov] = 0.0;
  prim[v3_cov] = 0.0;
  prim[WLORENTZ] = 1.0;

  if ((c2p.retain_B_atmo) && !(con[B1_con] != con[B1_con]) && !(con[B2_con] != con[B2_con])
        && !(con[B3_con] != con[B3_con])) {
    prim[B1_con] = con[B1_con];
    prim[B2_con] = con[B2_con];
    prim[B3_con] = con[B3_con];
  } else {
    prim[B1_con] = 0.0;
    prim[B2_con] = 0.0;
    prim[B3_con] = 0.0;
  }

  int keyerr = 0;
  int anyerr = 0;
  double xrho = prim[RHO];
  double xeps = 0.0;
  double xpress = 0.0;
  int keytemp;
  double xtemp,xye;

  if (c2p.evolve_Ye){
    prim[TEMP]  = c2p.T_atmo;
//    prim[TEMP]  = 0.1;
    prim[YE]    = c2p.Ye_atmo;
    prim[PRESS] = c2p.P_atmo;
    prim[EPS]   = c2p.eps_atmo;

    /* GB: I commented out this */
    /* // compute eps, press */
    /* keytemp = 1; */
    /* xtemp = prim[TEMP]; */
    /* xye = prim[YE]; */
    /* double xent = 0.0; */
    /* double xabar = 0.0; */

    /* EOS_press_ent_abar(c2p.eoskey,keytemp,xrho,&xeps,&xtemp,xye,&xpress,&xent,&xabar,&keyerr); */

    /* prim[ENT] = xent; */
    /* prim[A_BAR] = xabar; */

  } else {
    // compute eps, press
    keytemp = 0;
    xtemp = 0.0;
    xye = 0.0;

    EOS_Omni_press(c2p.eoskey_polytrope, keytemp, c2p.eos_prec, 1, &xrho, &xeps, &xtemp, &xye, &xpress, &keyerr, &anyerr);

    EOS_Omni_EpsFromPress(c2p.eoskey_polytrope, keytemp, c2p.eos_prec, 1, &xrho, &xeps, &xtemp, &xye, &xpress, &xeps, &keyerr, &anyerr);

    prim[EPS] = xeps;
    prim[PRESS] = xpress;

  }


}


void reset_velocity(const double gup[4][4], double * prim)
{

  /* Adjust initial guess for 3-velocity to ensure v^2 <= 1. Since the
   * metric is different from the previous time step, but the velocity components
   * as initial guesses are taken from the previous time step, v^2 <= 1 is not
   * guaranteed anymore with the present metric components.
   */

  double velupx = gup[1][1]*prim[v1_cov] + gup[1][2]*prim[v2_cov] + gup[1][3]*prim[v3_cov];
  double velupy = gup[2][1]*prim[v1_cov] + gup[2][2]*prim[v2_cov] + gup[2][3]*prim[v3_cov];
  double velupz = gup[3][1]*prim[v1_cov] + gup[3][2]*prim[v2_cov] + gup[3][3]*prim[v3_cov];

  double v2 = velupx*prim[v1_cov] + velupy*prim[v2_cov] + velupz*prim[v3_cov];

  if (v2 >= 1.0){
    prim[v1_cov] *= 1.0 / (v2*(1.0e0+1.0e-10));
    prim[v2_cov] *= 1.0 / (v2*(1.0e0+1.0e-10));
    prim[v3_cov] *= 1.0 / (v2*(1.0e0+1.0e-10));

    velupx = gup[1][1]*prim[v1_cov] + gup[1][2]*prim[v2_cov] + gup[1][3]*prim[v3_cov];
    velupy = gup[2][1]*prim[v1_cov] + gup[2][2]*prim[v2_cov] + gup[2][3]*prim[v3_cov];
    velupz = gup[3][1]*prim[v1_cov] + gup[3][2]*prim[v2_cov] + gup[3][3]*prim[v3_cov];

    v2 = velupx*prim[v1_cov] + velupy*prim[v2_cov] + velupz*prim[v3_cov];
    prim[WLORENTZ] = 1.0 / sqrt(1.0 - v2);
  } else if (v2 < 0.0){
    prim[v1_cov] = 0.0;
    prim[v2_cov] = 0.0;
    prim[v3_cov] = 0.0;
    prim[WLORENTZ] = 1.0;

  }
}

void reset_temperature(double * prim)
{
  // reset temperature to atmosphere level and update
  // affected primitives for consistency

  int keyerr = 0;
  int anyerr = 0;
  int keytemp = 1;
  double xrho = prim[RHO];
  double xeps = 0.0;
  double xpress = 0.0;
  double xtemp = c2p.T_atmo;
  double xye = prim[YE];


  if (c2p.evolve_T){
    double xent = 0.0;
    double xabar = 0.0;
    EOS_press_ent_abar(c2p.eoskey,keytemp,xrho,&xeps,&xtemp,xye,&xpress,&xent,&xabar,&keyerr);
    prim[ENT] = xent;
    prim[A_BAR] = xabar;
//    prim[MU_HAT] = 0.0;
  } else {
    EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &xrho, &xeps, &xtemp, &xye, &xpress, &keyerr, &anyerr);
  }

  // update primitives
  prim[EPS]   = xeps;
  prim[TEMP]  = xtemp;
  prim[PRESS] = xpress;
}


void c2p_recovery_scheme(const int c2p_method, struct c2p_report * c2p_rep,
      const double S_squared, const double BdotS, const double B_squared,
      const double * con, double * prim, const double original_prim[c2p.numprims],
      const double g_up[4][4], const double g_lo[4][4], const double tol_x, const int index, double xco, double yco,double zco,const int iteration, double * METRIC, double * METRIC_PHYS, double * METRIC_LAP_PSI4)
{

  // This routine selects the actual con2prim recovery scheme, any scheme can be called

  /*
  Note concerning Newton-Raphson root-finding schemes:
    Handling of initial guesses for primitives in con2prim
    SAFEGUESS = 0: use primitives from previous time step (default)
    SAFEGUESS = 1: use "safe guess" values from Cerdá-Durán 2008, Eq.(39)-(42)
    only used if SAFEGUESS=0 did not converge (see below)

  Note: stepsize parameter only needed for nuc EOS routines and Newton-Raphson root
    finding (although just a dummy parameter right now)

  */

  int retval; 
  switch(c2p_method) {
      case 1 :
        // 3D Newton-Raphson based on Cerdá-Durán 2008, Eq.(21),(22),(29) [method=1], i.e.
        // adding eps-eps(rho,T) = 0 as a third equation, or
        {
          int stepsize=1;
          int SAFEGUESS=0;
          NR_3D(c2p_rep,S_squared,BdotS,B_squared,con,prim,g_up,g_lo,
            tol_x,SAFEGUESS,stepsize);

          // restrict eps
          //if( (c2p_rep->failed) && (c2p.eoskey != 1) ){
          //  use_epsmin = true;
          //  // reset primitives (initial guesses)
          //  for(int l=0;l<c2p.numprims;l++) {
          //    prim[l] = original_prim[l];
          //  }
          //  NR_3D(c2p_rep,S_squared,BdotS,B_squared,con,prim,g_up,g_lo,
          //      tol_x,use_epsmin,SAFEGUESS,stepsize);
          //  use_epsmin = false;
          //}

          // retry with safe guess initial values
          // if(c2p_rep->failed){
          //   SAFEGUESS=1;
          //   use_epsmin = (c2p.eoskey == 1);
          //   NR_3D(c2p_rep,S_squared,BdotS,B_squared,con,prim,g_up,g_lo,
          //           tol_x,use_epsmin,SAFEGUESS,stepsize);
          //   if ((c2p_rep->failed) && (c2p.c2p_method == 1) && (c2p.eoskey != 1)){
          //     // try with restricted eps
          //    use_epsmin = true;
          //    NR_3D(c2p_rep,S_squared,BdotS,B_squared,con,prim,g_up,g_lo,
          //       tol_x,use_epsmin,SAFEGUESS,stepsize);
          //   }
          // }
          break;
        }
      case 2 :
        // 3D Newton-Raphson based on Cerdá-Durán 2008, Eq.(21),(22),(26) [method=2], i.e.
        // adding P-P(rho,T) = 0 as a third equation
        {
          int stepsize=1;
          int SAFEGUESS=0;
          NR_3D(c2p_rep,S_squared,BdotS,B_squared,con,prim,g_up,g_lo,
            tol_x,SAFEGUESS,stepsize);

          // retry with safe guess initial values
          if(c2p_rep->failed){
            SAFEGUESS=1;
            NR_3D(c2p_rep,S_squared,BdotS,B_squared,con,prim,g_up,g_lo,
              tol_x,SAFEGUESS,stepsize);
          }
          break;
        }
      case 3 :
        // variant of 2D Newton-Raphson based on Cerdá-Durán 2008, Eq.(21),(22) using
        // rest-mass density and temperature as root finding variables
        {
          int stepsize=1;
          int SAFEGUESS=0;
          NR_2D_rhoT(c2p_rep,S_squared,BdotS,B_squared,con,prim,g_lo,
              tol_x,SAFEGUESS,stepsize);

          // retry with safe guess initial values
          if(c2p_rep->failed){
            SAFEGUESS=1;
            NR_2D_rhoT(c2p_rep,S_squared,BdotS,B_squared,con,prim,g_lo,
              tol_x,SAFEGUESS,stepsize);
          }
          break;
        }
      case 4 :
        // variant of 2D Newton-Raphson based on Cerdá-Durán 2008, Eq.(21),(22) using
        // Lorentz factor and temperature as root finding variables
        {
          int stepsize=1;
          int SAFEGUESS=0;
          NR_2D_WT(c2p_rep,S_squared,BdotS,B_squared,con,prim,g_lo,
              tol_x,SAFEGUESS,stepsize);

          // retry with safe guess initial values
          if(c2p_rep->failed){
            SAFEGUESS=1;
            NR_2D_WT(c2p_rep,S_squared,BdotS,B_squared,con,prim,g_lo,
              tol_x,SAFEGUESS,stepsize);
          }
          break;
        }
      case 5 :
        {
          bool use_epsmin = 0;
          palenzuela(c2p_rep,S_squared,BdotS,B_squared,con,prim,g_up,g_lo,
            tol_x,use_epsmin,index,xco,yco,zco,iteration);
          break;
        }
      case 6 :
        {
          NR_2D_Noble(c2p_rep,S_squared,BdotS,B_squared,con,prim,g_up,g_lo,
            tol_x);
          break;
        }
      case 7 :
        {
          newman(c2p_rep,S_squared,BdotS,B_squared,con,prim,g_up,g_lo,tol_x,xco,yco,zco,iteration, METRIC_LAP_PSI4);
          break;
        }
      case 11 :
        {
          HARM(c2p_rep,S_squared,BdotS,B_squared,con,prim,g_up,g_lo,
            tol_x,xco,yco,zco,iteration, index, METRIC, METRIC_PHYS, METRIC_LAP_PSI4);
          break;
        }
      default :
        {
          printf("%s\n", "ERROR: No valid MHD con2prim scheme identifier!");
          exit(1);
          break;
        }
    }
}


void con2prim_MHD_(double * prim, double const * con, CCTK_REAL g_up[4][4],
      CCTK_REAL g_lo[4][4], const bool key_excise, const bool c2p_grace,
      struct c2p_report * c2p_rep, const int iteration, const int index, double xco, double yco, double zco, double * METRIC, double * METRIC_PHYS, double * METRIC_LAP_PSI4)
{

  bool reset2atmo = false;
  c2p_rep->failed = false;
  c2p_rep->retry  = false;
  c2p_rep->nEOScalls = 0;
  strcpy(c2p_rep->err_msg, "Everything all right!");

  double c2p_grace_fac = 1.0;
  double c2p_grace_atmo_fac = 1.0;
  if (c2p_grace) {
    c2p_grace_fac = 1e2;
    c2p_grace_atmo_fac = 4.0;
  }

  // check for NaNs in conservatives
  for(int l=0;l<c2p.numcons;l++){
    if (con[l] != con[l]){
      reset2atmo = true;
      c2p_rep->failed = true;
      strcpy( c2p_rep->err_msg, "NaN(s) in conservatives. Resetting to atmosphere.");
    }
  }

  // calculate some quantities for c2p scheme
  double B_squared, S_squared, BdotS, St_squared;
  calc_B2_S2_BdotS(con, g_up, g_lo, &B_squared, &S_squared, &St_squared, &BdotS, METRIC, METRIC_PHYS, METRIC_LAP_PSI4);

  if ((B_squared != B_squared) || (S_squared != S_squared) || (BdotS != BdotS))
  {
    reset2atmo = true;
    c2p_rep->failed = true;
    strcpy( c2p_rep->err_msg, "NaN(s) in B2, S2, or BdotS. Resetting to atmosphere.");
  }


  if ((reset2atmo) || (key_excise) || (con[D] <= c2p.rho_atmo * (1.0+c2p.rho_atmo_tol*c2p_grace_atmo_fac)))
  {
    // reset to atmosphere:
    // - if below density threshold
    // - if in excised region
    // - if NaNs in conservatives
    reset_to_atmosphere(prim, con);
    c2p_rep->adjust_cons = true;
    c2p_rep->count = 0;

  } else {

    // perform recovery of primitives

    // Adjust initial guess for 3-velocity to ensure v^2 <= 1.
    reset_velocity(g_up,prim);

    // keep original initial guesses to initialize retry (see below)
    double original_prim[c2p.numprims];
    for(int l=0;l<c2p.numprims;l++){
	    original_prim[l] = prim[l];
    }

    c2p_recovery_scheme(c2p.c2p_method,c2p_rep,S_squared,BdotS,B_squared,con,prim,original_prim,
		    g_up,g_lo,c2p.tol_x*c2p_grace_fac,index,xco,yco,zco,iteration,METRIC, METRIC_PHYS,METRIC_LAP_PSI4);
  
      if(c2p_rep->failed){

      // try backup scheme
      if (c2p.use_c2p_method_backup) {
        // reset primitives (initial guesses)
        for(int l=0;l<c2p.numprims;l++){
          prim[l] = original_prim[l];
        }
        c2p_recovery_scheme(c2p.c2p_method_backup,c2p_rep,S_squared,BdotS,B_squared,con,prim,original_prim,
          g_up,g_lo,c2p.tol_x*c2p_grace_fac,index,xco,yco,zco,iteration,METRIC, METRIC_PHYS,METRIC_LAP_PSI4);
      }

      // retry with larger tolerance
// #if !STANDALONE
//       if(c2p_rep->failed){
//         c2p_rep->retry = true;

//         // reset primitives (initial guesses)
//         for(int l=0;l<c2p.numprims;l++){
// 	        prim[l] = original_prim[l];
//         }

//         c2p_recovery_scheme(c2p.c2p_method,c2p_rep,S_squared,BdotS,B_squared,con,prim,original_prim,
//            g_up,g_lo,c2p.tol_x_retry*c2p_grace_fac);
//       }
// #endif
      if(c2p_rep->failed){
        // as a last resort, perform a prim2con and check whether we are actually
        // reasonably consistent with the conservatives anyway and accept
        // c2p attempt in this case
        //double con_check[c2p.numcons];
        //prim2con_MHD_(prim, con_check, g_up, g_lo);

        //c2p_rep->failed = false;
        //double deviation = 0.0;
        //for(int l=0;l<=7;l++){
        //  if (con[l] == 0.0){
        //    deviation = fabs(con_check[l]);
        //  }else{
        //    deviation = fabs((con_check[l] - con[l])/con[l]);
        //  }
        //  if ((deviation > 1e-5*c2p_grace_fac) || (deviation != deviation)){
        //    // con2prim really failed, reset to atmosphere
//            reset2atmo = true;
            c2p_rep->failed = true;
	    strcpy( c2p_rep->err_msg, "Retry failed. Trying Font fix.");
        //  }
        //}
      }
    }

    /*============================ START NEW ============================*/
    //before resetting to atmosphere, try a Font fix -- the Font fix routine
    //in the original Illinois solver has been generalized to in principle
    //allow for any EOS (checked for G=2 polytrope)
  
    // If it failed after the retry, use the new (generalized) Font fix subroutine 
    int font_fix_applied=0;
    if(c2p_rep->failed) {
        font_fix_applied=1;
        CCTK_REAL u_xl=1e100, u_yl=1e100, u_zl=1e100; // Set to insane values to ensure they are overwritten.
          int check = font_fix_general_EOS(&u_xl,&u_yl,&u_zl,con,prim,METRIC_PHYS,METRIC_LAP_PSI4);
  	c2p_rep->failed=false;
        //Translate to HARM primitive now:
	//these are already cavariant, so no need to change anything. In Illinois the Font fix returns covariant vector
	//and calculates a contravariant one from it which it puts into the prims array... no need to do that here
        prim[v1_cov] = u_xl;
        prim[v2_cov] = u_yl;
        prim[v3_cov] = u_zl;
        if (check==1) {
          CCTK_VInfo(CCTK_THORNSTRING,"Font fix failed!");
  	reset2atmo = true;
  	c2p_rep->failed=true;
        }
      }
    /*============================ END NEW ============================*/
    if ((reset2atmo) || (prim[RHO] <= c2p.rho_atmo * (1.0+c2p.rho_atmo_tol*c2p_grace_atmo_fac)))
    {
      // reset to atmosphere
      reset_to_atmosphere(prim, con);
      c2p_rep->adjust_cons = true;
      c2p_rep->failed=false;
    }

    if ((c2p.evolve_T) && (prim[TEMP] < c2p.T_atmo)) {
      // reset temperature if below temperature threshold
      reset_temperature(prim);
      c2p_rep->adjust_cons = true;
    }

    if (c2p.enforce_v2){
      // enforce that round off errors don't result in v2 >= 1!
      reset_velocity(g_up, prim);
    }
    /*============================ START NEW ============================*/
    //the calculation of post-recovery PRIMS is different between the old
    //Illnois version and these new solvers. Especially in the v-recovery
    //Illinois calculates the contravariant velocities and limits them
    //C2P calculates the covariant velocities, changes them to contravariant,
    //and doesn't limit them
    if(!c2p_rep->failed) {
      //Now that we have found some solution, we first limit velocity:
      //FIXME: Probably want to use exactly the same velocity limiter function here as in mhdflux.C
      //uti_new are contravariant velocities, so raise index using 3-metric
      CCTK_REAL utx_new = METRIC_PHYS[C2P_GUPXX]*prim[v1_cov] + METRIC_PHYS[C2P_GUPXY]*prim[v2_cov] + METRIC_PHYS[C2P_GUPXZ]*prim[v3_cov];
      CCTK_REAL uty_new = METRIC_PHYS[C2P_GUPXY]*prim[v1_cov] + METRIC_PHYS[C2P_GUPYY]*prim[v2_cov] + METRIC_PHYS[C2P_GUPYZ]*prim[v3_cov];
      CCTK_REAL utz_new = METRIC_PHYS[C2P_GUPXZ]*prim[v1_cov] + METRIC_PHYS[C2P_GUPYZ]*prim[v2_cov] + METRIC_PHYS[C2P_GUPZZ]*prim[v3_cov];
  
      //Velocity limiter:
      CCTK_REAL gijuiuj = METRIC_PHYS[C2P_GXX]*SQR(utx_new ) +
        2.0*METRIC_PHYS[C2P_GXY]*utx_new*uty_new + 2.0*METRIC_PHYS[C2P_GXZ]*utx_new*utz_new +
        METRIC_PHYS[C2P_GYY]*SQR(uty_new) + 2.0*METRIC_PHYS[C2P_GYZ]*uty_new*utz_new +
        METRIC_PHYS[C2P_GZZ]*SQR(utz_new);
      CCTK_REAL au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
      CCTK_REAL u0L = (au0m1+1.0)*METRIC_LAP_PSI4[C2P_LAPSEINV];
      CCTK_REAL lap = METRIC_LAP_PSI4[C2P_LAPSE];
  
      CCTK_REAL GAMMA_SPEED_LIMIT=10.0;
      // *** Limit velocity
      if (au0m1 > 0.9999999*(GAMMA_SPEED_LIMIT-1.0)) {
        CCTK_REAL fac = sqrt((SQR(GAMMA_SPEED_LIMIT)-1.0)/(SQR(1.0+au0m1) - 1.0));
        utx_new *= fac;
        uty_new *= fac;
        utz_new *= fac;
        gijuiuj = gijuiuj * SQR(fac);
        au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
        // Reset rho_b and u0
        u0L = (au0m1+1.0)*METRIC_LAP_PSI4[C2P_LAPSEINV];
        prim[RHO] =  con[D]/(METRIC_LAP_PSI4[C2P_LAPSE]*u0L);
        prim[YECON] =  con[YE]/(METRIC_LAP_PSI4[C2P_LAPSE]*u0L);
      } //Finished limiting velocity
  
    //The Font fix only sets the velocities.  Here we set the pressure & density HARM primitives.
      if(font_fix_applied==1) {
        prim[RHO] = con[D]/(METRIC_LAP_PSI4[C2P_LAPSE]*u0L);
        prim[YECON] = con[YE]/(METRIC_LAP_PSI4[C2P_LAPSE]*u0L);
        //Next set P = P_cold:
        const int keytemp = 0;
	const int eoskey = 4;
        double xeps = prim[PRESS]/prim[RHO];
        double xye = con[YE]/con[D];
        double xrho = prim[RHO];
        double xtemp = 0.1;
	double P_cold = 0.0;
        int anyerr = 0;
        int keyerr = 0;
        EOS_press_cold(&xrho, &xye, &P_cold, &xeps, &keyerr, &anyerr);
        prim[PRESS] = P_cold;
        prim[EPS] = xeps;
      } //Finished setting remaining primitives if there was a Font fix.

      prim[v1_cov] = (METRIC_PHYS[C2P_GXX]*utx_new + METRIC_PHYS[C2P_GXY]*uty_new + METRIC_PHYS[C2P_GXZ]*utz_new)/(u0L*lap);
      prim[v2_cov] = (METRIC_PHYS[C2P_GXY]*utx_new + METRIC_PHYS[C2P_GYY]*uty_new + METRIC_PHYS[C2P_GYZ]*utz_new)/(u0L*lap);
      prim[v3_cov] = (METRIC_PHYS[C2P_GXZ]*utx_new + METRIC_PHYS[C2P_GYZ]*uty_new + METRIC_PHYS[C2P_GZZ]*utz_new)/(u0L*lap);
  
    /*============================ END NEW ============================*/
    }
  }
}

/*============================ START NEW ============================*/
//this is the same Font fix routine for rho and vi in Illinois, but generalized from a Gamma-law (assuming a COLD EOS)
inline int font_fix_general_EOS(CCTK_REAL *u_x, CCTK_REAL *u_y, CCTK_REAL *u_z,CCTK_REAL *con,CCTK_REAL *prim,CCTK_REAL *METRIC_PHYS,CCTK_REAL *METRIC_LAP_PSI4) {
  CCTK_REAL rhob;
  CCTK_REAL tol = 1.e-15;
  CCTK_REAL Bxbar = prim[B1_con]*ONE_OVER_SQRT_4PI;
  CCTK_REAL Bybar = prim[B2_con]*ONE_OVER_SQRT_4PI;
  CCTK_REAL Bzbar = prim[B3_con]*ONE_OVER_SQRT_4PI;
  CCTK_REAL Bbar_x = METRIC_PHYS[C2P_GXX]*Bxbar + METRIC_PHYS[C2P_GXY]*Bybar + METRIC_PHYS[C2P_GXZ]*Bzbar;
  CCTK_REAL Bbar_y = METRIC_PHYS[C2P_GXY]*Bxbar + METRIC_PHYS[C2P_GYY]*Bybar + METRIC_PHYS[C2P_GYZ]*Bzbar;
  CCTK_REAL Bbar_z = METRIC_PHYS[C2P_GXZ]*Bxbar + METRIC_PHYS[C2P_GYZ]*Bybar + METRIC_PHYS[C2P_GZZ]*Bzbar;
  CCTK_REAL B2bar = Bxbar*Bbar_x + Bybar*Bbar_y + Bzbar*Bbar_z;
  CCTK_REAL Bbar = sqrt(B2bar);
  CCTK_REAL check_B_small = fabs(Bxbar)+fabs(Bybar)+fabs(Bzbar);
  if (check_B_small>0 && check_B_small<1.e-150) {
    // need to compute B2bar specially to prevent floating-point underflow
    CCTK_REAL Bmax = fabs(Bxbar);
    if (Bmax < fabs(Bybar)) Bmax=fabs(Bybar);
    if (Bmax < fabs(Bzbar)) Bmax=fabs(Bzbar);
    CCTK_REAL Bxtmp=Bxbar/Bmax, Bytemp=Bybar/Bmax, Bztemp=Bzbar/Bmax;
    CCTK_REAL B_xtemp=Bbar_x/Bmax, B_ytemp=Bbar_y/Bmax, B_ztemp=Bbar_z/Bmax;
    Bbar = sqrt(Bxtmp*B_xtemp + Bytemp*B_ytemp + Bztemp*B_ztemp)*Bmax;
  }

  CCTK_REAL BbardotS = Bxbar*con[S1_cov]*METRIC_LAP_PSI4[C2P_PSI6] + Bybar*con[S2_cov]*METRIC_LAP_PSI4[C2P_PSI6] + Bzbar*con[S3_cov]*METRIC_LAP_PSI4[C2P_PSI6];
  CCTK_REAL BbardotS2 = BbardotS*BbardotS;
  CCTK_REAL hatBbardotS = BbardotS/Bbar;
  if (Bbar<1.e-300) hatBbardotS = 0.0;
  CCTK_REAL Psim6 = 1.0/METRIC_LAP_PSI4[C2P_PSI6];

  CCTK_REAL sdots = (METRIC_PHYS[C2P_GUPXX]*SQR(con[S1_cov]) + METRIC_PHYS[C2P_GUPYY]*SQR(con[S2_cov]) + METRIC_PHYS[C2P_GUPZZ]*SQR(con[S3_cov])
    + 2.0*( METRIC_PHYS[C2P_GUPXY]*con[S1_cov]*con[S2_cov] + METRIC_PHYS[C2P_GUPXZ]*con[S1_cov]*con[S3_cov]
            + METRIC_PHYS[C2P_GUPYZ]*con[S2_cov]*con[S3_cov]))*METRIC_LAP_PSI4[C2P_PSI6]*METRIC_LAP_PSI4[C2P_PSI6];

  if (sdots<1.e-300) {
    rhob = con[D]*METRIC_LAP_PSI4[C2P_PSI6]*Psim6;
    *u_x=0.0; *u_y=0.0; *u_z=0.0;
    return 0;
  }

  // Initial guess for W, S_fluid and rhob
  CCTK_REAL W0 = sqrt( SQR(hatBbardotS) + SQR(con[D]*METRIC_LAP_PSI4[C2P_PSI6]) ) * Psim6;
  CCTK_REAL Sf20 = (SQR(W0)*sdots + BbardotS2*(B2bar + 2.0*W0))/SQR(W0+B2bar);
  CCTK_REAL rhob0 = con[D]*METRIC_LAP_PSI4[C2P_PSI6]*Psim6/sqrt(1.0+Sf20/SQR(con[D]*METRIC_LAP_PSI4[C2P_PSI6]));
  CCTK_REAL W=W0,Sf2=Sf20,rhob1=rhob0;
  // Assume cold for Font fix?
  CCTK_REAL xye = con[YE]/con[D];
  CCTK_INT keytemp = 1;
  //ASSUMES COLD EOS, NO REASON NOT TO PASS IT THE RECOVERED TEMP
  CCTK_REAL xtemp = 0.1;

  //****************************************************************
  //                          FONT FIX
  // Impose Font fix when HARM primitives solver fails to find
  //   acceptable set of primitives.
  //****************************************************************
  bool fontcheck=true;

  int itcount = 0, maxits=500;
  while(fontcheck && itcount < maxits) {
    itcount++;
    W0 = W;
    Sf20 = Sf2;
    rhob0 = rhob1;

    // first, find rhob for the given S_fluid^2
    CCTK_REAL xeps; CCTK_REAL xpress; int anyerr = 0;int keyerr = 0;
    EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &rhob0, &xeps, &xtemp, &xye, &xpress, &keyerr, &anyerr);
    CCTK_REAL enth = 1.0 + xeps + xpress/rhob0;
    rhob1 = con[D]*METRIC_LAP_PSI4[C2P_PSI6]*Psim6/sqrt(1.0+Sf20/SQR(con[D]*METRIC_LAP_PSI4[C2P_PSI6]*enth));
    while( fabs(rhob1-rhob0) > rhob1*tol) {
      rhob0 = rhob1;
      //calculate enthalpy for a general EOS
      CCTK_REAL xeps; CCTK_REAL xpress; int anyerr = 0;int keyerr = 0;
      EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &rhob0, &xeps, &xtemp, &xye, &xpress, &keyerr, &anyerr);
      CCTK_REAL enth = 1.0 + xeps + xpress/rhob0;
      rhob1 = con[D]*METRIC_LAP_PSI4[C2P_PSI6]*Psim6/sqrt(1.0+Sf20/SQR(con[D]*METRIC_LAP_PSI4[C2P_PSI6]*enth));
    }

    //enth = 1.0 + xeps + xP/rhob0;
    anyerr = 0; keyerr = 0;
    EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &rhob0, &xeps, &xtemp, &xye, &xpress, &keyerr, &anyerr);
    enth = 1.0 + xeps + xpress/rhob0;
    W = sqrt( Sf20 + SQR(con[D]*METRIC_LAP_PSI4[C2P_PSI6]*enth))*Psim6;
    Sf2 = (SQR(W)*sdots + BbardotS2*(B2bar + 2.0*W))/SQR(W+B2bar);
    if ( fabs(W-W0) < W*tol && fabs(Sf20-Sf2) < Sf2*tol) fontcheck=false;
  }
  
  if (itcount>=maxits) {
    // Increase tol and try again 
    maxits*=100;
    tol *=10.0;
    itcount = 0;
    fontcheck=true;
    while(fontcheck && itcount < maxits) {
      itcount++;
      W0 = W;
      Sf20 = Sf2;
      rhob0 = rhob1;

      // first, find rhob for the given S_fluid^2
      CCTK_REAL xeps; CCTK_REAL xpress; int anyerr = 0;int keyerr = 0;
      EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &rhob0, &xeps, &xtemp, &xye, &xpress, &keyerr, &anyerr);
      CCTK_REAL enth = 1.0 + xeps + xpress/rhob0;
      rhob1 = con[D]*METRIC_LAP_PSI4[C2P_PSI6]*Psim6/sqrt(1.0+Sf20/SQR(con[D]*METRIC_LAP_PSI4[C2P_PSI6]*enth));
      while( fabs(rhob1-rhob0) > rhob1*tol) {
        rhob0 = rhob1;
        anyerr = 0; keyerr = 0;
        EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &rhob0, &xeps, &xtemp, &xye, &xpress, &keyerr, &anyerr);
        enth = 1.0 + xeps + xpress/rhob0;
        rhob1 = con[D]*METRIC_LAP_PSI4[C2P_PSI6]*Psim6/sqrt(1.0+Sf20/SQR(con[D]*METRIC_LAP_PSI4[C2P_PSI6]*enth));
      }

      anyerr = 0; keyerr = 0;
      EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &rhob0, &xeps, &xtemp, &xye, &xpress, &keyerr, &anyerr);
      enth = 1.0 + xeps + xpress/rhob0;
      W = sqrt( Sf20 + SQR(con[D]*METRIC_LAP_PSI4[C2P_PSI6]*enth))*Psim6;
      Sf2 = (SQR(W)*sdots + BbardotS2*(B2bar + 2.0*W))/SQR(W+B2bar);
      if ( fabs(W-W0) < W*tol && fabs(Sf20-Sf2) < Sf2*tol) fontcheck=false;
    }
  }
  //************************************************************************************************************** 

  if(fontcheck==true) {
    return 1;
  }

  // Font fix works, now compute u_i
  rhob = rhob1;
  CCTK_REAL xeps; CCTK_REAL xpress; int anyerr = 0;int keyerr = 0;
  EOS_Omni_press(c2p.eoskey, keytemp, c2p.eos_prec, 1, &rhob0, &xeps, &xtemp, &xye, &xpress, &keyerr, &anyerr);
  CCTK_REAL enth = 1.0 + xeps + xpress/rhob0;
  CCTK_REAL gammav = con[D]*METRIC_LAP_PSI4[C2P_PSI6]*Psim6/rhob;
  CCTK_REAL rhosh = con[D]*METRIC_LAP_PSI4[C2P_PSI6]*enth;
  CCTK_REAL fac1 = METRIC_LAP_PSI4[C2P_PSI6]*BbardotS/(gammav*rhosh);
  CCTK_REAL fac2 = 1.0/(rhosh + METRIC_LAP_PSI4[C2P_PSI6]*B2bar/gammav);
  *u_x = fac2*(con[S1_cov]*METRIC_LAP_PSI4[C2P_PSI6] + fac1*Bbar_x);
  *u_y = fac2*(con[S2_cov]*METRIC_LAP_PSI4[C2P_PSI6] + fac1*Bbar_y);
  *u_z = fac2*(con[S3_cov]*METRIC_LAP_PSI4[C2P_PSI6] + fac1*Bbar_z);
  return 0;
}
/*============================ END NEW ============================*/
