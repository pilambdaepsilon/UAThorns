//-----------------------------------------------------------------------------
// Compute the flux for advecting rho_star, tau (Font's energy variable),
//  and S_i .
//-----------------------------------------------------------------------------
static inline void mhdflux(const int eos_key,const CCTK_REAL c2pprec,int i,int j,int k,const int flux_dirn,CCTK_REAL *Ul,CCTK_REAL *Ur,  CCTK_REAL *FACEVAL,CCTK_REAL *FACEVAL_LAPSE_PSI4,eos_struct &eos,
                           CCTK_REAL &cmax,CCTK_REAL &cmin,
                           CCTK_REAL &rho_star_flux,CCTK_REAL &tau_flux,CCTK_REAL &st_x_flux,CCTK_REAL &st_y_flux,CCTK_REAL &st_z_flux,CCTK_REAL &Yet_flux) {


  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL psi4 = FACEVAL_LAPSE_PSI4[PSI4];
  CCTK_REAL psi6 = FACEVAL_LAPSE_PSI4[PSI4]*FACEVAL_LAPSE_PSI4[PSI2];
  CCTK_REAL psim4 = 1.0/(psi4);

  CCTK_REAL alpha_sqrt_gamma = FACEVAL_LAPSE_PSI4[LAPSE]*psi6;
  CCTK_REAL ONE_OVER_LAPSE = 1.0/FACEVAL_LAPSE_PSI4[LAPSE];
  CCTK_REAL ONE_OVER_LAPSE_SQUARED=SQR(ONE_OVER_LAPSE);

  // First compute P_{cold}, \epsilon_{cold}, dP_{cold}/drho, \epsilon_{th}, h, and \Gamma_{cold},
  // for right and left faces:
  // of these, we only need h_{l,r} for calculation of the fluxes (the others are used in calculation of cs2, which 
  // has been generalized to all EOSs)
  CCTK_REAL P_coldr=0.0,eps_coldr=0.0,dPcold_drhor=0,eps_thr=0,h_r=0,gamma_coldr=0.0;
  CCTK_REAL P_coldl=0.0,eps_coldl=0.0,dPcold_drhol=0,eps_thl=0,h_l=0,gamma_coldl=0.0;
  CCTK_REAL P_l=0.0; CCTK_REAL P_r=0.0;
  if(!use_ConservativeToPrimitive){
    compute_P_cold__eps_cold__dPcold_drho__eps_th__h__gamma_cold(Ur,eos,P_coldr,eps_coldr,dPcold_drhor,eps_thr,h_r,gamma_coldr);
    compute_P_cold__eps_cold__dPcold_drho__eps_th__h__gamma_cold(Ul,eos,P_coldl,eps_coldl,dPcold_drhol,eps_thl,h_l,gamma_coldl);
    P_l = Ul[PRESSURE]; P_r = Ur[PRESSURE];
  }
  else{
    if(eos_key==4){
       /******************* ENF ******************/
       //For consistency, use the reconstructed eps, rho, and Ye to get
       //P and T, then use these to get the enthalpy. currently
       //there is no function to get eps from pressure in the EOS driver
       //left face
       CCTK_REAL xrhol   = Ul[RHOB]; CCTK_REAL xepsl   = Ul[EPS]; CCTK_REAL xyel    = Ul[YE];
       CCTK_REAL xtempl  = Ul[TEMP]; CCTK_REAL xpressl  = Ul[PRESSURE];    // just initializing
       CCTK_INT  keytemp = 0;     // We enter with eps (reconstructed)
       CCTK_INT  npoints  = 1; CCTK_INT  anyerr, keyerr;
       EOS_Omni_press(eos_key,keytemp,c2pprec,npoints,&xrhol,&xepsl,&xtempl,&xyel,&xpressl,&keyerr,&anyerr);
       //right face
       CCTK_REAL xrhor   = Ur[RHOB]; CCTK_REAL xepsr   = Ur[EPS]; CCTK_REAL xyer    = Ur[YE];
       CCTK_REAL xtempr  = Ur[TEMP]; CCTK_REAL xpressr  = Ur[PRESSURE];    // just initializing
       keytemp = 0;     // We enter with eps (reconstructed)
       npoints  = 1;
       EOS_Omni_press(eos_key,keytemp,c2pprec,npoints,&xrhor,&xepsr,&xtempr,&xyer,&xpressr,&keyerr,&anyerr);
       P_l = xpressl;    P_r = xpressr;
       h_r = 1.0 + xepsr + P_r/xrhor;
       h_l = 1.0 + xepsl + P_l/xrhol;
  
       /******************* ETHFIX -- NOT NEEDED BUT CAN HELP ******************/
       // eps_th fix: if eps < eps_cold, append a thermal gamma-law component to the pressure
       // as this suggests that eps_th < 0.0 which is nonsense for positive temperatures
       // first, get the cold pressure and epsilon (r and l)
       CCTK_REAL P_coldr=0.0,eps_coldr=0.0,eps_thr=0.0;
       CCTK_REAL P_coldl=0.0,eps_coldl=0.0,eps_thl=0;
  
       EOS_press_cold(&xrhol, &xyel, &P_coldl, &eps_coldl, &keyerr, &anyerr);
       EOS_press_cold(&xrhor, &xyer, &P_coldr, &eps_coldr, &keyerr, &anyerr);
       // then, get the 'thermal' epsilon, and append a 'thermal' pressure if needed
       eps_thl = xepsl - eps_coldl;
       eps_thr = xepsr - eps_coldr;
       if(eps_thl < 0.0){P_l = P_coldl + (gamma_th-1)*xrhol*eps_thl;h_l = 1.0 + xepsl + P_l/xrhol;}
       if(eps_thr < 0.0){P_r = P_coldr + (gamma_th-1)*xrhor*eps_thr;h_r = 1.0 + xepsr + P_r/xrhor;}
    }
    else{
       /******************* REC ******************/
       // this is what is typically done (in compute_P_cold__eps_cold__dPcold_drho__eps_th__h__gamma_cold), 
       // but it leads to inconsistencies for general (finite temperature) EOSs, as reconstructed EPS 
       // need not be consistent with reconstructed PRESSURE. This is not troublesome for hybrid EOSs 
       // and fixes the ocassional crash caused by not setting eps=eps_cold in 
       // apply_tau_floor__enforce_limits_on_primitives_and_recompute_conservs.C
       
       P_l=Ul[PRESSURE]; P_r=Ur[PRESSURE];
       h_r = 1.0 + Ur[EPS] + Ur[PRESSURE]/Ur[RHOB];
       h_l = 1.0 + Ul[EPS] + Ul[PRESSURE]/Ul[RHOB];
    }
  }

  //Compute face velocities
  // Begin by computing u0
  output_stats stats; stats.failure_checker=0;
  CCTK_REAL u0_r,u0_l;
  impose_speed_limit_output_u0(FACEVAL,Ur,psi4,ONE_OVER_LAPSE,stats,u0_r);
  impose_speed_limit_output_u0(FACEVAL,Ul,psi4,ONE_OVER_LAPSE,stats,u0_l);

  //Next compute b^{\mu}, the magnetic field measured in the comoving fluid frame:
  CCTK_REAL ONE_OVER_LAPSE_SQRT_4PI = ONE_OVER_LAPSE*ONE_OVER_SQRT_4PI;
  /***********************************************************/
  /********** RIGHT FACE ************/
  // Note that smallbr[4] = b^a defined in Gammie's paper, on the right face.
  CCTK_REAL u_x_over_u0_psi4r,u_y_over_u0_psi4r,u_z_over_u0_psi4r;
  CCTK_REAL smallbr[NUMVARS_SMALLB];
  // Compute b^{a}, b^2, and u_i over u^0
  compute_smallba_b2_and_u_i_over_u0_psi4(FACEVAL,FACEVAL_LAPSE_PSI4,Ur,u0_r,ONE_OVER_LAPSE_SQRT_4PI,
                                          u_x_over_u0_psi4r,u_y_over_u0_psi4r,u_z_over_u0_psi4r,smallbr);
  // Then compute u_xr,u_yr, and u_zr. We need to set the zeroth component so we can specify U_LOWER{r,l}[{UX,UY,UZ}] (UX=1,UY=2,UZ=3).
  CCTK_REAL U_LOWERr[4] = { 0.0, u_x_over_u0_psi4r*u0_r*FACEVAL_LAPSE_PSI4[PSI4], u_y_over_u0_psi4r*u0_r*FACEVAL_LAPSE_PSI4[PSI4],
                            u_z_over_u0_psi4r*u0_r*FACEVAL_LAPSE_PSI4[PSI4] };
  /********** LEFT FACE ************/
  // Note that smallbl[4] = b^a defined in Gammie's paper, on the left face.
  CCTK_REAL u_x_over_u0_psi4l,u_y_over_u0_psi4l,u_z_over_u0_psi4l;
  CCTK_REAL smallbl[NUMVARS_SMALLB];
  // Compute b^{a}, b^2, and u_i over u^0
  compute_smallba_b2_and_u_i_over_u0_psi4(FACEVAL,FACEVAL_LAPSE_PSI4,Ul,u0_l,ONE_OVER_LAPSE_SQRT_4PI,
                                          u_x_over_u0_psi4l,u_y_over_u0_psi4l,u_z_over_u0_psi4l,smallbl);
  // Then compute u_xr,u_yr, and u_zr. We need to set the zeroth component so we can specify U_LOWER{r,l}[{UX,UY,UZ}]
  CCTK_REAL U_LOWERl[4] = { 0.0, u_x_over_u0_psi4l*u0_l*FACEVAL_LAPSE_PSI4[PSI4], u_y_over_u0_psi4l*u0_l*FACEVAL_LAPSE_PSI4[PSI4],
                            u_z_over_u0_psi4l*u0_l*FACEVAL_LAPSE_PSI4[PSI4] };
  /***********************************************************/

  // Compute v02 = v_A^2 + c_s^2*(1.0-v_A^2), where c_s = sound speed, and v_A = Alfven velocity
  CCTK_REAL v02r,v02l;
  // First right face
  compute_v02(eos_key,c2pprec,dPcold_drhor,eos.gamma_th,eps_thr,h_r,smallbr,Ur,v02r);
  // Then left face.
  compute_v02(eos_key,c2pprec,dPcold_drhol,eos.gamma_th,eps_thl,h_l,smallbl,Ul,v02l);

  int offset=flux_dirn-1;

  // Now that we have computed v02 = v_A^2 + c_s^2*(1.0-v_A^2), we can
  //   next compute c_+ and c_- using a simplified dispersion relation.
  //   Note that, in solving the simplified disp. relation, we overestimate
  //   c_+ and c_- by around a factor of 2, making the MHD evolution more
  //   diffusive (and potentially more *stable*) than it could be.
  CCTK_REAL cplusr,cminusr,cplusl,cminusl;
  find_cp_cm(cplusr,cminusr,v02r,u0_r,
             Ur[VX+offset],ONE_OVER_LAPSE_SQUARED,FACEVAL[SHIFTX+offset],psim4,FACEVAL[GUPXX+offset]);
  find_cp_cm(cplusl,cminusl,v02l,u0_l,
             Ul[VX+offset],ONE_OVER_LAPSE_SQUARED,FACEVAL[SHIFTX+offset],psim4,FACEVAL[GUPXX+offset]);

  // Then compute cmax, cmin. This is required for the HLL flux.
  CCTK_REAL cmaxL =  MAX(0.0,MAX(cplusl,cplusr));
  CCTK_REAL cminL = -MIN(0.0,MIN(cminusl,cminusr));

  //*********************************************************************
  // density flux = \rho_* v^m, where m is the current flux direction (the m index)
  //*********************************************************************
  CCTK_REAL rho_star_r = alpha_sqrt_gamma*Ur[RHOB]*u0_r;
  CCTK_REAL rho_star_l = alpha_sqrt_gamma*Ul[RHOB]*u0_l;
  CCTK_REAL Yet_r = rho_star_r*Ur[YE];
  CCTK_REAL Yet_l = rho_star_l*Ul[YE];
  CCTK_REAL Fr = rho_star_r*Ur[VX+offset]; // flux_dirn = 2, so offset = 1, implies Ur[VX] -> Ur[VY]
  CCTK_REAL Fl = rho_star_l*Ul[VX+offset]; // flux_dirn = 2, so offset = 1, implies Ul[VX] -> Ul[VY]
  CCTK_REAL Fr_Yet = Yet_r*Ur[VX+offset];
  CCTK_REAL Fl_Yet = Yet_l*Ul[VX+offset];

  // HLL step for rho_star:
  rho_star_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(rho_star_r-rho_star_l) )/(cmaxL + cminL);
  Yet_flux = (cminL*Fr_Yet + cmaxL*Fl_Yet - cminL*cmaxL*(Yet_r-Yet_l) )/(cmaxL + cminL);

  //*********************************************************************
  // energy flux = \alpha^2 \sqrt{\gamma} T^{0m} - \rho_* v^m, where m is the current flux direction (the m index)
  //*********************************************************************
  // First compute some useful metric quantities:
  CCTK_REAL alpha_squared_sqrt_gamma = FACEVAL_LAPSE_PSI4[LAPSE]*alpha_sqrt_gamma;
  CCTK_REAL g4uptm =  ONE_OVER_LAPSE_SQUARED*FACEVAL[SHIFTX+offset];
  CCTK_REAL g4uptt = -ONE_OVER_LAPSE_SQUARED;
  /********** RIGHT FACE ************/
  // Compute a couple useful hydro quantities:
  CCTK_REAL rho0_h_plus_b2_r = (Ur[RHOB]*h_r + smallbr[SMALLB2]);
  CCTK_REAL P_plus_half_b2_r = (P_r+0.5*smallbr[SMALLB2]);
  // Then compute T^{0m} and the flux:
  CCTK_REAL TUP0m_r = rho0_h_plus_b2_r*SQR(u0_r)*Ur[VX+offset] + P_plus_half_b2_r*g4uptm - smallbr[SMALLBT]*smallbr[SMALLBX+offset];
  Fr = alpha_squared_sqrt_gamma * TUP0m_r - rho_star_r * Ur[VX+offset];
  // Finally compute tau
  CCTK_REAL TUP00_r = rho0_h_plus_b2_r*u0_r*u0_r + P_plus_half_b2_r*g4uptt - smallbr[SMALLBT]*smallbr[SMALLBT];
  CCTK_REAL tau_r = alpha_squared_sqrt_gamma * TUP00_r - rho_star_r;
  /********** LEFT FACE *************/
  // Compute a couple useful hydro quantities:
  CCTK_REAL rho0_h_plus_b2_l = (Ul[RHOB]*h_l + smallbl[SMALLB2]);
  CCTK_REAL P_plus_half_b2_l = (P_l+0.5*smallbl[SMALLB2]);
  // Then compute T^{0m} and the flux:
  CCTK_REAL TUP0m_l = rho0_h_plus_b2_l*SQR(u0_l)*Ul[VX+offset] + P_plus_half_b2_l*g4uptm - smallbl[SMALLBT]*smallbl[SMALLBX+offset];
  Fl = alpha_squared_sqrt_gamma * TUP0m_l - rho_star_l * Ul[VX+offset];
  // Finally compute tau
  CCTK_REAL TUP00_l = rho0_h_plus_b2_l*u0_l*u0_l + P_plus_half_b2_l*g4uptt - smallbl[SMALLBT]*smallbl[SMALLBT];
  CCTK_REAL tau_l = alpha_squared_sqrt_gamma * TUP00_l - rho_star_l;

  // HLL step for tau:
  tau_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(tau_r-tau_l) )/(cmaxL + cminL);

  //*********************************************************************
  // momentum flux = \alpha \sqrt{\gamma} T^m_j, where m is the current flux direction (the m index)
  //*********************************************************************
  // b_j = g_{ij} (b^i + b^t shift^i), g_{ij} = physical metric
  //CCTK_REAL sbtr=0,sbtl=0;
  CCTK_REAL smallb_lowerr[NUMVARS_SMALLB],smallb_lowerl[NUMVARS_SMALLB];
  lower_4vector_output_spatial_part(psi4,FACEVAL,smallbr,smallb_lowerr);
  lower_4vector_output_spatial_part(psi4,FACEVAL,smallbl,smallb_lowerl);

  /********** Flux for S_x **********/
  // [S_x flux] = \alpha \sqrt{\gamma} T^m_x, where m is the current flux direction (the m index)
  //    Again, offset = 0 for reconstruction in x direction, 1 for y, and 2 for z
  //    Note that kronecker_delta[flux_dirn][0] = { 1 if flux_dirn==1, 0 otherwise }.
  Fr = alpha_sqrt_gamma*( rho0_h_plus_b2_r*(u0_r*Ur[VX+offset])*U_LOWERr[UX]
                          + P_plus_half_b2_r*kronecker_delta[flux_dirn][0] - smallbr[SMALLBX+offset]*smallb_lowerr[SMALLBX] );
  Fl = alpha_sqrt_gamma*( rho0_h_plus_b2_l*(u0_l*Ul[VX+offset])*U_LOWERl[UX]
                          + P_plus_half_b2_l*kronecker_delta[flux_dirn][0] - smallbl[SMALLBX+offset]*smallb_lowerl[SMALLBX] );

  //        S_x =\alpha\sqrt{\gamma}( T^0_x )
  CCTK_REAL st_x_r = alpha_sqrt_gamma*( rho0_h_plus_b2_r*u0_r*U_LOWERr[UX] - smallbr[SMALLBT]*smallb_lowerr[SMALLBX] );
  CCTK_REAL st_x_l = alpha_sqrt_gamma*( rho0_h_plus_b2_l*u0_l*U_LOWERl[UX] - smallbl[SMALLBT]*smallb_lowerl[SMALLBX] );

  // HLL step for Sx:
  st_x_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(st_x_r-st_x_l) )/(cmaxL + cminL);

  /********** Flux for S_y **********/
  // [S_y flux] = \alpha \sqrt{\gamma} T^m_y, where m is the current flux direction (the m index)
  //    Again, offset = 1 for reconstruction in x direction, 2 for y, and 3 for z
  //    Note that kronecker_delta[flux_dirn][1] = { 1 if flux_dirn==2, 0 otherwise }.
  Fr = alpha_sqrt_gamma*( rho0_h_plus_b2_r*(u0_r*Ur[VX+offset])*U_LOWERr[UY] + P_plus_half_b2_r*kronecker_delta[flux_dirn][1]
                          - smallbr[SMALLBX+offset]*smallb_lowerr[SMALLBY] );
  Fl = alpha_sqrt_gamma*( rho0_h_plus_b2_l*(u0_l*Ul[VX+offset])*U_LOWERl[UY] + P_plus_half_b2_l*kronecker_delta[flux_dirn][1]
                          - smallbl[SMALLBX+offset]*smallb_lowerl[SMALLBY] );

  //        S_y =\alpha\sqrt{\gamma}( T^0_y )
  CCTK_REAL st_y_r = alpha_sqrt_gamma*( rho0_h_plus_b2_r*u0_r*U_LOWERr[UY] - smallbr[SMALLBT]*smallb_lowerr[SMALLBY] );
  CCTK_REAL st_y_l = alpha_sqrt_gamma*( rho0_h_plus_b2_l*u0_l*U_LOWERl[UY] - smallbl[SMALLBT]*smallb_lowerl[SMALLBY] );

  // HLL step for Sy:
  st_y_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(st_y_r-st_y_l) )/(cmaxL + cminL);
  /********** Flux for S_z **********/
  // [S_z flux] = \alpha \sqrt{\gamma} T^m_z, where m is the current flux direction (the m index)
  //    Again, offset = 1 for reconstruction in x direction, 2 for y, and 3 for z
  //    Note that kronecker_delta[flux_dirn][2] = { 1 if flux_dirn==3, 0 otherwise }.
  Fr = alpha_sqrt_gamma*( rho0_h_plus_b2_r*(u0_r*Ur[VX+offset])*U_LOWERr[UZ] + P_plus_half_b2_r*kronecker_delta[flux_dirn][2]
                          - smallbr[SMALLBX+offset]*smallb_lowerr[SMALLBZ] );
  Fl = alpha_sqrt_gamma*( rho0_h_plus_b2_l*(u0_l*Ul[VX+offset])*U_LOWERl[UZ] + P_plus_half_b2_l*kronecker_delta[flux_dirn][2]
                          - smallbl[SMALLBX+offset]*smallb_lowerl[SMALLBZ] );

  //        S_z =\alpha\sqrt{\gamma}( T^0_z )
  CCTK_REAL st_z_r = alpha_sqrt_gamma*( rho0_h_plus_b2_r*u0_r*U_LOWERr[UZ] - smallbr[SMALLBT]*smallb_lowerr[SMALLBZ] );
  CCTK_REAL st_z_l = alpha_sqrt_gamma*( rho0_h_plus_b2_l*u0_l*U_LOWERl[UZ] - smallbl[SMALLBT]*smallb_lowerl[SMALLBZ] );

  // HLL step for Sz:
  st_z_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(st_z_r-st_z_l) )/(cmaxL + cminL);

  cmax = cmaxL;
  cmin = cminL;
}
