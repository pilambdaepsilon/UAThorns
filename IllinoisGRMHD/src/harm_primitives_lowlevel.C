inline int harm_primitives_gammalaw_lowlevel(const int index,const int i,const int j,const int k,CCTK_REAL *X,CCTK_REAL *Y,CCTK_REAL *Z,
                                             CCTK_REAL *METRIC,CCTK_REAL *METRIC_PHYS,CCTK_REAL *METRIC_LAP_PSI4,
                                             CCTK_REAL *CONSERVS,CCTK_REAL *PRIMS,
                                             CCTK_REAL g4dn[NDIM][NDIM],CCTK_REAL g4up[NDIM][NDIM],
                                             struct output_stats &stats, eos_struct &eos) {

  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL kpoly = eos.k_tab[0];
  CCTK_REAL gamma = eos.gamma_tab[0];
  int gamma_equals2 = 1;
  if (fabs(gamma-2.0) > 1.e-10) gamma_equals2 = 0;

  // declare some variables for HARM.
  CCTK_REAL U[NPR];						//array with (#of primitives) indices to hold conservatives 
  CCTK_REAL prim[NPR];						// " primitives
  CCTK_REAL detg = METRIC_LAP_PSI4[LAPSE]*METRIC_LAP_PSI4[PSI6]; // == alpha sqrt{gamma} = alpha Psi^6

  // Check to see if the metric is positive-definite.
  // Note that this will slow down the code, and if the metric doesn't obey this, the run is probably too far gone to save,
  //   though if it happens deep in the horizon, it might resurrect the run.
  /*
    CCTK_REAL lam1,lam2,lam3;
    CCTK_REAL M11 = METRIC[GXX], M12=METRIC[GXY], M13=METRIC[GXZ], M22=METRIC[GYY], M23=METRIC[GYZ], M33=METRIC[GZZ];
    eigenvalues_3by3_real_sym_matrix(lam1, lam2, lam3,M11, M12, M13, M22, M23, M33);
    if (lam1 < 0.0 || lam2 < 0.0 || lam3 < 0.0) {
    // Metric is not positive-defitive, reset the physical metric to be conformally-flat.
    METRIC_PHYS[GXX] = METRIC_LAP_PSI4[PSI4];
    METRIC_PHYS[GXY] = 0.0;
    METRIC_PHYS[GXZ] = 0.0;
    METRIC_PHYS[GYY] = METRIC_LAP_PSI4[PSI4];
    METRIC_PHYS[GYZ] = 0.0;
    METRIC_PHYS[GZZ] = METRIC_LAP_PSI4[PSI4];
    METRIC_PHYS[GUPXX] = METRIC_LAP_PSI4[PSIM4];
    METRIC_PHYS[GUPXY] = 0.0;
    METRIC_PHYS[GUPXZ] = 0.0;
    METRIC_PHYS[GUPYY] = METRIC_LAP_PSI4[PSIM4];
    METRIC_PHYS[GUPYZ] = 0.0;
    METRIC_PHYS[GUPZZ] = METRIC_LAP_PSI4[PSIM4];
    }
  */

  // Note that ONE_OVER_SQRT_4PI gets us to the object
  // referred to as B^i in the Noble et al paper (and
  // apparently also in the comments to their code).
  // This is NOT the \mathcal{B}^i, which differs by 
  // a factor of the lapse.
  CCTK_REAL BxL_over_alpha_sqrt_fourpi = PRIMS[BX_CENTER]*METRIC_LAP_PSI4[LAPSEINV]*ONE_OVER_SQRT_4PI;
  CCTK_REAL ByL_over_alpha_sqrt_fourpi = PRIMS[BY_CENTER]*METRIC_LAP_PSI4[LAPSEINV]*ONE_OVER_SQRT_4PI;
  CCTK_REAL BzL_over_alpha_sqrt_fourpi = PRIMS[BZ_CENTER]*METRIC_LAP_PSI4[LAPSEINV]*ONE_OVER_SQRT_4PI;

  CCTK_REAL rho_b_oldL = PRIMS[RHOB];
  CCTK_REAL YeL        = PRIMS[YE];
  CCTK_REAL P_oldL     = PRIMS[PRESSURE];
  CCTK_REAL vxL        = PRIMS[VX];
  CCTK_REAL vyL        = PRIMS[VY];
  CCTK_REAL vzL        = PRIMS[VZ];

  /*
    -- Driver for new prim. var. solver.  The driver just translates
    between the two sets of definitions for U and P.  The user may 
    wish to alter the translation as they see fit.  


    //         /  rho u^t           \                             //
    //    U =  |  T^t_t   + rho u^t |  sqrt(-det(g_{\mu\nu}))     //
    //         |  T^t_i             |                             //
    //         \   B^i              /                             //
    //                                                            //
    //        /    rho        \                                   //
    //    P = |    uu         |                                   //
    //        | \tilde{u}^i   |                                   //
    //        \   B^i         /                                   //

    (above equations have been fixed by Yuk Tung & Zach)
  */
  
  // U[NPR]    = conserved variables (current values on input/output);
  // g4dn[NDIM][NDIM] = covariant form of the 4-metric ;
  // g4up[NDIM][NDIM] = contravariant form of the 4-metric ;
  // gdet             = sqrt( - determinant of the 4-metric) ;
  // prim[NPR] = primitive variables (guess on input, calculated values on 
  //                     output if there are no problems);

  // U[1]   =  
  // U[2-4] =  stildei + rhostar

  CCTK_REAL rho_star_orig = CONSERVS[RHOSTAR];
  CCTK_REAL mhd_st_x_orig = CONSERVS[STILDEX];
  CCTK_REAL mhd_st_y_orig = CONSERVS[STILDEY];
  CCTK_REAL mhd_st_z_orig = CONSERVS[STILDEZ];
  CCTK_REAL tau_orig      = CONSERVS[TAUENERGY];
  CCTK_REAL Yet_orig      = CONSERVS[YET];

  // Other ideas for setting the gamma speed limit
  //CCTK_REAL GAMMA_SPEED_LIMIT = 100.0;
  //if(METRIC_LAP_PSI4[PSI6]>Psi6threshold) GAMMA_SPEED_LIMIT=500.0;
  //if(METRIC_LAP_PSI4[PSI6]>Psi6threshold) GAMMA_SPEED_LIMIT=100.0;

  //FIXME: Only works if poisoning is turned on. Otherwise will access unknown memory. This trick alone speeds up the whole code (Cowling) by 2%.
  //int startguess=0;
  //if(std::isnan(PRIMS[VX])) startguess=1;
  int startguess=1;

  CCTK_REAL u0L=1.0;

  for(int which_guess=startguess;which_guess<3;which_guess++) {
    int check;

    if(which_guess==1) {
      //Use a different initial guess:
      rho_b_oldL = CONSERVS[RHOSTAR]/METRIC_LAP_PSI4[PSI6];
      YeL = CONSERVS[YET]/METRIC_LAP_PSI4[PSI6];
      if (gamma_equals2==1) {
        P_oldL = kpoly*rho_b_oldL*rho_b_oldL;
      } else {
        P_oldL = kpoly*pow(rho_b_oldL,gamma);
      }
      u0L = METRIC_LAP_PSI4[LAPSEINV];
      vxL = -METRIC[SHIFTX];
      vyL = -METRIC[SHIFTY];
      vzL = -METRIC[SHIFTZ];
    }

    if(which_guess==2) {
      //Use atmosphere as initial guess:
      rho_b_oldL = 100.0*rho_b_atm;
//      YeL = 100.0*rho_b_atm;
//      YeL = 0.5;
      // GAMMA=2 ONLY:
      if (gamma_equals2==1) {
        P_oldL = kpoly*rho_b_oldL*rho_b_oldL;
      } else {
        P_oldL = kpoly*pow(rho_b_oldL,gamma);
      }
      u0L = METRIC_LAP_PSI4[LAPSEINV];
      vxL = -METRIC[SHIFTX];
      vyL = -METRIC[SHIFTY];
      vzL = -METRIC[SHIFTZ];
    }

    // Fill the array of conserved variables according to the wishes of Utoprim_2d.
    U[RHO]    = CONSERVS[RHOSTAR];
    U[YECON]    = CONSERVS[YET];
    U[UU]     = -CONSERVS[TAUENERGY]*METRIC_LAP_PSI4[LAPSE] - (METRIC_LAP_PSI4[LAPSE]-1.0)*CONSERVS[RHOSTAR] + 
      METRIC[SHIFTX]*CONSERVS[STILDEX] + METRIC[SHIFTY]*CONSERVS[STILDEY]  + METRIC[SHIFTZ]*CONSERVS[STILDEZ] ; // note the minus sign on tau
    U[UTCON1] = CONSERVS[STILDEX];
    U[UTCON2] = CONSERVS[STILDEY];
    U[UTCON3] = CONSERVS[STILDEZ];
    U[BCON1]  = detg*BxL_over_alpha_sqrt_fourpi;
    U[BCON2]  = detg*ByL_over_alpha_sqrt_fourpi;
    U[BCON3]  = detg*BzL_over_alpha_sqrt_fourpi;

    CCTK_REAL uL = P_oldL/(gamma_th /* <- Should be local polytropic gamma factor */ - 1.0);
    CCTK_REAL utxL = u0L*(vxL + METRIC[SHIFTX]);
    CCTK_REAL utyL = u0L*(vyL + METRIC[SHIFTY]);
    CCTK_REAL utzL = u0L*(vzL + METRIC[SHIFTZ]);

    prim[RHO]    = rho_b_oldL;
    prim[YECON]    = YeL;
    prim[UU]     = uL;
    prim[UTCON1] = utxL;
    prim[UTCON2] = utyL;
    prim[UTCON3] = utzL;
    prim[BCON1]  = BxL_over_alpha_sqrt_fourpi;
    prim[BCON2]  = ByL_over_alpha_sqrt_fourpi;
    prim[BCON3]  = BzL_over_alpha_sqrt_fourpi;

    /*************************************************************/
    // CALL HARM PRIMITIVES SOLVER:
    check = Utoprim_2d(U, g4dn, g4up, detg, prim,stats.n_iter);
    // Note that we have modified this solver, so that nearly 100% 
    // of the time it yields either a good root, or a root with 
    // negative epsilon (i.e., pressure).
    /*************************************************************/

    // Use the new Font fix subroutine 
    int font_fix_applied=0;
    if(check!=0) {
      font_fix_applied=1;
      CCTK_REAL u_xl=1e100, u_yl=1e100, u_zl=1e100; // Set to insane values to ensure they are overwritten.
      if (gamma_equals2==1) {
        check = font_fix_gamma_equals2(u_xl,u_yl,u_zl,CONSERVS,PRIMS,METRIC_PHYS,METRIC_LAP_PSI4,eos);
      } else {
        check = font_fix_general_gamma(u_xl,u_yl,u_zl,CONSERVS,PRIMS,METRIC_PHYS,METRIC_LAP_PSI4,eos);
      }
      //Translate to HARM primitive now:
      prim[UTCON1] = METRIC_PHYS[GUPXX]*u_xl + METRIC_PHYS[GUPXY]*u_yl + METRIC_PHYS[GUPXZ]*u_zl;
      prim[UTCON2] = METRIC_PHYS[GUPXY]*u_xl + METRIC_PHYS[GUPYY]*u_yl + METRIC_PHYS[GUPYZ]*u_zl;
      prim[UTCON3] = METRIC_PHYS[GUPXZ]*u_xl + METRIC_PHYS[GUPYZ]*u_yl + METRIC_PHYS[GUPZZ]*u_zl;
      if (check==1) {
        CCTK_VInfo(CCTK_THORNSTRING,"Font fix failed!");
        CCTK_VInfo(CCTK_THORNSTRING,"i,j,k = %d %d %d, stats.failure_checker = %d x,y,z = %e %e %e , index=%d st_i = %e %e %e, rhostar = %e, Yet = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",i,j,k,stats.failure_checker,X[index],Y[index],Z[index],index,mhd_st_x_orig,mhd_st_y_orig,mhd_st_z_orig,rho_star_orig,Yet_orig,PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER],METRIC_PHYS[GXX],METRIC_PHYS[GXY],METRIC_PHYS[GXZ],METRIC_PHYS[GYY],METRIC_PHYS[GYZ],METRIC_PHYS[GZZ],METRIC_LAP_PSI4[PSI6]);
      }
    }
    stats.failure_checker+=font_fix_applied*10000;
    stats.font_fixed=font_fix_applied;
    /*************************************************************/

    if(check==0) {
      //Now that we have found some solution, we first limit velocity:
      //FIXME: Probably want to use exactly the same velocity limiter function here as in mhdflux.C
      CCTK_REAL utx_new = prim[UTCON1];
      CCTK_REAL uty_new = prim[UTCON2];
      CCTK_REAL utz_new = prim[UTCON3];

      //Velocity limiter:
      CCTK_REAL gijuiuj = METRIC_PHYS[GXX]*SQR(utx_new ) +
        2.0*METRIC_PHYS[GXY]*utx_new*uty_new + 2.0*METRIC_PHYS[GXZ]*utx_new*utz_new +
        METRIC_PHYS[GYY]*SQR(uty_new) + 2.0*METRIC_PHYS[GYZ]*uty_new*utz_new +
        METRIC_PHYS[GZZ]*SQR(utz_new);
      CCTK_REAL au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
      u0L = (au0m1+1.0)*METRIC_LAP_PSI4[LAPSEINV];

      // *** Limit velocity
      stats.vel_limited=0;
      if (au0m1 > 0.9999999*(GAMMA_SPEED_LIMIT-1.0)) {
        CCTK_REAL fac = sqrt((SQR(GAMMA_SPEED_LIMIT)-1.0)/(SQR(1.0+au0m1) - 1.0));
        utx_new *= fac;
        uty_new *= fac;
        utz_new *= fac;
        gijuiuj = gijuiuj * SQR(fac);
        au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
        // Reset rho_b and u0
        u0L = (au0m1+1.0)*METRIC_LAP_PSI4[LAPSEINV];
        prim[RHO] =  rho_star_orig/(METRIC_LAP_PSI4[LAPSE]*u0L*METRIC_LAP_PSI4[PSI6]);
        prim[YECON] =  Yet_orig/(METRIC_LAP_PSI4[LAPSE]*u0L*METRIC_LAP_PSI4[PSI6]);
        stats.vel_limited=1;
        stats.failure_checker+=1000;
      } //Finished limiting velocity

      //The Font fix only sets the velocities.  Here we set the pressure & density HARM primitives.
      if(font_fix_applied==1) {
        prim[RHO] = rho_star_orig/(METRIC_LAP_PSI4[LAPSE]*u0L*METRIC_LAP_PSI4[PSI6]);
        prim[YECON] = Yet_orig/(METRIC_LAP_PSI4[LAPSE]*u0L*METRIC_LAP_PSI4[PSI6]);
        //Next set P = P_cold:
        CCTK_REAL P_cold;
        if (gamma_equals2==1) {
          P_cold = kpoly*prim[RHO]*prim[RHO]; // <- FIXME: GAMMA=2 ONLY
        } else {
          P_cold = kpoly*pow(prim[RHO],gamma);
        }
        prim[UU] = P_cold/(gamma_th /* <- Should be local polytropic gamma factor */-1.0);
      } //Finished setting remaining primitives if there was a Font fix.

      PRIMS[RHOB]     = prim[RHO];
      PRIMS[YE]       = Yet_orig/rho_star_orig;
      PRIMS[PRESSURE] = (gamma_th /* <- Should be local polytropic gamma factor */-1.0)*prim[UU];
      //Already set u0L.
      PRIMS[VX]       = utx_new/u0L - METRIC[SHIFTX];
      PRIMS[VY]       = uty_new/u0L - METRIC[SHIFTY];
      PRIMS[VZ]       = utz_new/u0L - METRIC[SHIFTZ];

      return 0;
    } else {
      //If we didn't find a root, then try again with a different guess.
    }
  }
  CCTK_VInfo(CCTK_THORNSTRING,"Couldn't find root from: [tau]%e [rho_star]%e [Yet]%e [st_i]%e %e %e, rhob approx=%e, rho_b_atm=%e, Bx=%e, By=%e, Bz=%e, gij_phys=%e %e %e %e %e %e, alpha=%e",
	     tau_orig,rho_star_orig,Yet_orig,mhd_st_x_orig,mhd_st_y_orig,mhd_st_z_orig,rho_star_orig/METRIC_LAP_PSI4[PSI6],rho_b_atm,PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER],METRIC_PHYS[GXX],METRIC_PHYS[GXY],METRIC_PHYS[GXZ],METRIC_PHYS[GYY],METRIC_PHYS[GYZ],METRIC_PHYS[GZZ],METRIC_LAP_PSI4[LAPSE]);
  return 1;
}

#include "harm_u2p_util.c"
#include "harm_utoprim_2d.c"
#include "eigen.C"
#include "font_fix_gamma_law.C"
