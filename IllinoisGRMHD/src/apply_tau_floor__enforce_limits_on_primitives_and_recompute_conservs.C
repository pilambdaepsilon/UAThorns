void eigenvalues_3by3_real_sym_matrix(CCTK_REAL & lam1, CCTK_REAL & lam2, CCTK_REAL & lam3,
                                      CCTK_REAL M11, CCTK_REAL M12, CCTK_REAL M13, CCTK_REAL M22, CCTK_REAL M23, CCTK_REAL M33);

static inline int apply_tau_floor(const int index,const CCTK_REAL tau_atm,const CCTK_REAL rho_b_atm,const CCTK_REAL Psi6threshold,CCTK_REAL *PRIMS,CCTK_REAL *METRIC,CCTK_REAL *METRIC_PHYS,CCTK_REAL *METRIC_LAP_PSI4,output_stats &stats,eos_struct &eos,  CCTK_REAL *CONSERVS) {

  CCTK_REAL gamma = eos.gamma_tab[0];
  CCTK_REAL kpoly = eos.k_tab[0];
  int gamma_equals2 = 1;
  if (fabs(gamma-2.0) > 1.e-10) gamma_equals2 = 0;

  //First apply the rho_star floor:

  //rho_star = alpha u0 Psi6 rho_b, alpha u0 > 1, so if rho_star < Psi6 rho_b_atm, then we are GUARANTEED that we can reset to atmosphere.
  //if(CONSERVS[RHOSTAR] < 1e4*METRIC_LAP_PSI4[PSI6]*rho_b_atm) {
  //if(CONSERVS[RHOSTAR] < 2*METRIC_LAP_PSI4[PSI6]*rho_b_atm) {

  CCTK_REAL lam1,lam2,lam3;
  eigenvalues_3by3_real_sym_matrix(lam1, lam2, lam3,METRIC[GXX], METRIC[GXY], METRIC[GXZ], METRIC[GYY], METRIC[GYZ], METRIC[GZZ]);
  if (lam1 < 0.0 || lam2 < 0.0 || lam3 < 0.0) {
    // Metric is not positive-defitive, reset the metric to be conformally-flat.
    METRIC[GXX] = 1.0;
    METRIC[GXY] = 0.0;
    METRIC[GXZ] = 0.0;
    METRIC[GYY] = 1.0;
    METRIC[GYZ] = 0.0;
    METRIC[GZZ] = 1.0;
    METRIC_PHYS[GUPXX] = METRIC_LAP_PSI4[PSIM4];
    METRIC_PHYS[GUPXY] = 0.0;
    METRIC_PHYS[GUPXZ] = 0.0;
    METRIC_PHYS[GUPYY] = METRIC_LAP_PSI4[PSIM4];
    METRIC_PHYS[GUPYZ] = 0.0;
    METRIC_PHYS[GUPZZ] = METRIC_LAP_PSI4[PSIM4];
  }

  //Next, prepare for the tau and stilde fixes:

  CCTK_REAL Bxbar = PRIMS[BX_CENTER]*ONE_OVER_SQRT_4PI,Bybar = PRIMS[BY_CENTER]*ONE_OVER_SQRT_4PI,Bzbar = PRIMS[BZ_CENTER]*ONE_OVER_SQRT_4PI;

  CCTK_REAL Bbar_x = METRIC_PHYS[GXX]*Bxbar + METRIC_PHYS[GXY]*Bybar + METRIC_PHYS[GXZ]*Bzbar;
  CCTK_REAL Bbar_y = METRIC_PHYS[GXY]*Bxbar + METRIC_PHYS[GYY]*Bybar + METRIC_PHYS[GYZ]*Bzbar;
  CCTK_REAL Bbar_z = METRIC_PHYS[GXZ]*Bxbar + METRIC_PHYS[GYZ]*Bybar + METRIC_PHYS[GZZ]*Bzbar;
  CCTK_REAL Bbar2 = Bxbar*Bbar_x + Bybar*Bbar_y + Bzbar*Bbar_z;
  CCTK_REAL Bbar = sqrt(Bbar2);
  CCTK_REAL check_B_small = fabs(Bxbar)+fabs(Bybar)+fabs(Bzbar);
  if (check_B_small>0 && check_B_small<1.e-150) {
    // need to compute Bbar specially to prevent floating-point underflow
    CCTK_REAL Bmax = fabs(Bxbar);
    if (Bmax < fabs(Bybar)) Bmax=fabs(Bybar);
    if (Bmax < fabs(Bzbar)) Bmax=fabs(Bzbar);
    CCTK_REAL Bxtmp=Bxbar/Bmax, Bytemp=Bybar/Bmax, Bztemp=Bzbar/Bmax;
    CCTK_REAL B_xtemp=Bbar_x/Bmax, B_ytemp=Bbar_y/Bmax, B_ztemp=Bbar_z/Bmax;
    Bbar = sqrt(Bxtmp*B_xtemp + Bytemp*B_ytemp + Bztemp*B_ztemp)*Bmax;
  }
  CCTK_REAL BbardotS = Bxbar*CONSERVS[STILDEX] + Bybar*CONSERVS[STILDEY] + Bzbar*CONSERVS[STILDEZ];
  CCTK_REAL hatBbardotS = BbardotS/Bbar;
  if (Bbar<1.e-300) hatBbardotS = 0.0;

  // Limit hatBbardotS
  //CCTK_REAL max_gammav = 100.0;
  //CCTK_REAL rhob_max = CONSERVS[RHOSTAR]/METRIC_LAP_PSI4[PSI6];
  //CCTK_REAL hmax = 1.0 + 2.0*rhob_max;
  //CCTK_REAL abs_hatBbardotS_max = sqrt(SQR(max_gammav)-1.0)*CONSERVS[RHOSTAR]*hmax;
  //if (fabs(hatBbardotS) > abs_hatBbardotS_max) {
  //   CCTK_REAL fac_reduce = abs_hatBbardotS_max/fabs(hatBbardotS);
  //   CCTK_REAL hatBbardotS_max = hatBbardotS*fac_reduce;
  //   CCTK_REAL Bbar_inv = 1.0/Bbar;
  //   CCTK_REAL hat_Bbar_x = Bbar_x*Bbar_inv;
  //   CCTK_REAL hat_Bbar_y = Bbar_y*Bbar_inv;
  //   CCTK_REAL hat_Bbar_z = Bbar_z*Bbar_inv;
  //   CCTK_REAL sub_fact = hatBbardotS_max - hatBbardotS;
  //   CONSERVS[STILDEX] += sub_fact*hat_Bbar_x;
  //   CONSERVS[STILDEY] += sub_fact*hat_Bbar_y;
  //   CONSERVS[STILDEZ] += sub_fact*hat_Bbar_z;
  //   hatBbardotS = hatBbardotS_max;
  //   BbardotS *= fac_reduce;
  //   CONSERVS[STILDEX] = CONSERVS[STILDEX]; CONSERVS[STILDEY] = CONSERVS[STILDEY]; CONSERVS[STILDEZ] = CONSERVS[STILDEZ];
  //}

  CCTK_REAL sdots= METRIC_PHYS[GUPXX]*SQR(CONSERVS[STILDEX])+METRIC_PHYS[GUPYY]*SQR(CONSERVS[STILDEY])+METRIC_PHYS[GUPZZ]*SQR(CONSERVS[STILDEZ])+2.0*
    (METRIC_PHYS[GUPXY]*CONSERVS[STILDEX]*CONSERVS[STILDEY]+METRIC_PHYS[GUPXZ]*CONSERVS[STILDEX]*CONSERVS[STILDEZ]+METRIC_PHYS[GUPYZ]*CONSERVS[STILDEY]*CONSERVS[STILDEZ]);

  CCTK_REAL Wm = sqrt(SQR(hatBbardotS)+ SQR(CONSERVS[RHOSTAR]))/METRIC_LAP_PSI4[PSI6];
  CCTK_REAL Sm2 = (SQR(Wm)*sdots + SQR(BbardotS)*(Bbar2+2.0*Wm))/SQR(Wm+Bbar2);
  CCTK_REAL Wmin = sqrt(Sm2 + SQR(CONSERVS[RHOSTAR]))/METRIC_LAP_PSI4[PSI6];
  CCTK_REAL sdots_fluid_max = sdots;

  //tau fix, applicable when B==0 and B!=0:
  if(CONSERVS[TAUENERGY] < 0.5*METRIC_LAP_PSI4[PSI6]*Bbar2) {
    CONSERVS[TAUENERGY] = tau_atm+0.5*METRIC_LAP_PSI4[PSI6]*Bbar2;
    stats.failure_checker+=1000000;
  }

  CCTK_REAL tau_fluid_min = CONSERVS[TAUENERGY] - 0.5*METRIC_LAP_PSI4[PSI6]*Bbar2 - (Bbar2*sdots - SQR(BbardotS))*0.5/(METRIC_LAP_PSI4[PSI6]*SQR(Wmin+Bbar2));

  //Apply Stilde fix when B==0.
  //if(PRIMS[BX_CENTER]==0 && PRIMS[BY_CENTER]==0 && PRIMS[BZ_CENTER]==0 && (METRIC_LAP_PSI4[PSI6]>30.0 || CONSERVS[RHOSTAR]/METRIC_LAP_PSI4[PSI6]<100*rho_b_atm)) {
  //if(check_B_small < 1.e-300) {
  CCTK_REAL Patm;
  if (gamma_equals2==1) {
    Patm = kpoly*rho_b_atm*rho_b_atm;
  } else {
    Patm = kpoly*pow(rho_b_atm,gamma);
  }
  if(check_B_small*check_B_small < Patm*1e-32) {
    CCTK_REAL rhot=CONSERVS[TAUENERGY]*(CONSERVS[TAUENERGY]+2.0*CONSERVS[RHOSTAR]);
    CCTK_REAL safetyfactor = 0.999999;
    //if(METRIC_LAP_PSI4[PSI6]>Psi6threshold) safetyfactor=0.99;

    if(sdots > safetyfactor*rhot) {
      CCTK_REAL rfactm1 = sqrt((safetyfactor*rhot)/sdots);
      CONSERVS[STILDEX]*=rfactm1;
      CONSERVS[STILDEY]*=rfactm1;
      CONSERVS[STILDEZ]*=rfactm1;
      stats.failure_checker+=10000000;
    }
  } else if(METRIC_LAP_PSI4[PSI6]>Psi6threshold) {
    //Apply new Stilde fix.
    if (tau_fluid_min < tau_atm*1.001) {
      tau_fluid_min = tau_atm*1.001;
      CONSERVS[TAUENERGY] = tau_fluid_min + 0.5*METRIC_LAP_PSI4[PSI6]*Bbar2 + (Bbar2*sdots - SQR(BbardotS))*0.5/(METRIC_LAP_PSI4[PSI6]*SQR(Wmin+Bbar2));
    }
    CCTK_REAL LHS = tau_fluid_min*(tau_fluid_min+2.0*CONSERVS[RHOSTAR]);
    CCTK_REAL RHS = sdots_fluid_max;

    CCTK_REAL safetyfactor = 0.999999;
    if(safetyfactor*LHS < RHS) {
      CCTK_REAL rfactm1 = sqrt((safetyfactor*LHS)/RHS);
      CONSERVS[STILDEX]*=rfactm1;
      CONSERVS[STILDEY]*=rfactm1;
      CONSERVS[STILDEZ]*=rfactm1;
      stats.failure_checker+=100000000;
    }
  }



  return 0;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/


void IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(const int already_computed_physical_metric_and_inverse,CCTK_REAL *U,struct output_stats &stats,eos_struct &eos,
                                                                       CCTK_REAL *METRIC,CCTK_REAL g4dn[4][4],CCTK_REAL g4up[4][4], CCTK_REAL *TUPMUNU,CCTK_REAL *TDNMUNU,CCTK_REAL *CONSERVS, const int index) {
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL METRIC_LAP_PSI4[NUMVARS_METRIC_AUX];
  SET_LAPSE_PSI4(METRIC_LAP_PSI4,METRIC);

  // Useful debugging tool, can be used to track fixes:
  //CCTK_REAL rho_b_orig=U[RHOB],P_orig=U[PRESSURE],vx_orig=U[VX],vy_orig=U[VY],vz_orig=U[VZ];

  /***********************************************************/
  // Enforce limits on pressure, density, and v^i
  /***********************************************************/
  // Density floor:
  U[RHOB] = MAX(U[RHOB],rho_b_atm);
  // Density ceiling:
  U[RHOB] = MIN(U[RHOB],rho_b_max);
//  U[EPS] = MAX(U[EPS],rho_b_atm);
//  U[EPS] = MIN(U[EPS],rho_b_max);
  // Ye floor and ceiling:
  U[YE] = MAX(U[YE],0.036);
  U[YE] = MIN(U[YE],1.0);
  U[TEMP] = MAX(U[TEMP],0.1);
  //Entropy floor
  U[ENTROPY] = MAX(U[ENTROPY],rho_b_atm);





  //   Next set h, the enthalpy:
  CCTK_REAL h_enthalpy=0,  P_cold=0,eps_cold=0,dPcold_drho,eps_th,gamma_cold; /* <- Note that in setting h, we need to define several
                                                                               *    other variables. Though some will be unused later
                                                                               *    in this function, they may be useful in other
                                                                               *    functions */
  CCTK_REAL eps;

  if (!use_ConservativeToPrimitive){
    compute_P_cold__eps_cold__dPcold_drho__eps_th__h__gamma_cold(U,eos,P_cold,eps_cold,dPcold_drho,eps_th,h_enthalpy,gamma_cold);
    // Pressure floor & ceiling:
    enforce_pressure_floor_ceiling(stats,eos.K_poly,P_cold,METRIC_LAP_PSI4[PSI6],Psi6threshold,U[RHOB],rho_b_atm,  U[PRESSURE]);
    // Possibly adjusted pressure, so recompute eps & h:
    eps = eps_cold + (U[PRESSURE]-P_cold)/(gamma_th-1.0)/U[RHOB];
  }else{
    eps = U[EPS];
    CCTK_REAL xrho = U[RHOB];
    CCTK_REAL xye  = U[YE];
    CCTK_REAL xeps;
    CCTK_INT anyerr, keyerr;
    CCTK_INT npoints = 1;
    // I have to consider the smallest temperature available to get the "cold" component
    EOS_press_cold(&xrho, &xye, &P_cold, &xeps, &keyerr, &anyerr);
    // Pressure floor & ceiling:
    enforce_pressure_floor_ceiling(stats,eos.K_poly,P_cold,METRIC_LAP_PSI4[PSI6],Psi6threshold,U[RHOB],rho_b_atm, U[PRESSURE]);
    // FIXME: I am not recomputing the correct epsilon here!
    // Hopefully floor and ceilings are applied only in interior of black holes and in the atmosphere
    // NOTE: Setting eps=xeps (returned from EOS_press_cold) here forces a cold enthalpy, which leads 
    // to incorrect evolution for finite Temp EOS. However, not setting eps=xeps sometimes leads to a crash 
    // for hybrid EOS. 
    // There is no function to get eps from pressure (without the temperature) in the 
    // EOS driver, so we leave eps = U[EPS] (recovered), even if it may be inconsistent with the pressure.
    // This is fixed by using the reconstructed pressure and epsilon to get the enthalpy in mhdflux.C
    //eps=xeps;
  }

  h_enthalpy = 1.0 + eps + U[PRESSURE]/U[RHOB];

  CCTK_REAL uUP[4];
  impose_speed_limit_output_u0(METRIC,U,METRIC_LAP_PSI4[PSI4],METRIC_LAP_PSI4[LAPSEINV],stats, uUP[0]);
  // Compute u^i. We've already set uUP[0] in the lines above.
  for(int ii=0;ii<3;ii++) uUP[UX+ii] = uUP[0]*U[VX+ii];

  // Useful debugging tool, part 2: can be used to track fixes:
  //if(P_orig!=U[PRESSURE] || rho_b_orig!=U[RHOB] || vx_orig!=U[VX] || vy_orig!=U[VY] || vz_orig!=U[VZ]) {

  /***************************************************************/
  // COMPUTE TUPMUNU, TDNMUNU, AND  CONSERVATIVES FROM PRIMITIVES
  /***************************************************************/
  // Compute b^{\mu}, b^2, and u_i/(u^0 Psi4)
  CCTK_REAL ONE_OVER_LAPSE_SQRT_4PI = METRIC_LAP_PSI4[LAPSEINV]*ONE_OVER_SQRT_4PI;
  CCTK_REAL u_x_over_u0_psi4,u_y_over_u0_psi4,u_z_over_u0_psi4;
  CCTK_REAL smallb[NUMVARS_SMALLB];
  compute_smallba_b2_and_u_i_over_u0_psi4(METRIC,METRIC_LAP_PSI4,U,uUP[0],ONE_OVER_LAPSE_SQRT_4PI,
                                          u_x_over_u0_psi4,u_y_over_u0_psi4,u_z_over_u0_psi4,smallb);
  // Compute u_i; we compute u_0 below.
  CCTK_REAL uDN[4] = { 1e200, u_x_over_u0_psi4*uUP[0]*METRIC_LAP_PSI4[PSI4],u_y_over_u0_psi4*uUP[0]*METRIC_LAP_PSI4[PSI4],u_z_over_u0_psi4*uUP[0]*METRIC_LAP_PSI4[PSI4] };

  // Precompute some useful quantities, for later:
  CCTK_REAL alpha_sqrt_gamma=METRIC_LAP_PSI4[LAPSE]*METRIC_LAP_PSI4[PSI6];
  CCTK_REAL rho0_h_plus_b2 = (U[RHOB]*h_enthalpy + smallb[SMALLB2]);
  CCTK_REAL P_plus_half_b2 = (U[PRESSURE]+0.5*smallb[SMALLB2]);

  if(!already_computed_physical_metric_and_inverse) {
    // If you like, see Eq 2.119 in Numerical Relativity, by Baumgarte & Shapiro:
    CCTK_REAL ONE_OVER_LAPSE_SQUARED = SQR(METRIC_LAP_PSI4[LAPSEINV]);

    // g^{\mu \nu} = upper four-metric.
    //CCTK_REAL g4up[4][4];
    g4up[0][0] = -ONE_OVER_LAPSE_SQUARED;
    g4up[0][1] = ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTX];
    g4up[0][2] = ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTY];
    g4up[0][3] = ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTZ];
    g4up[1][1] = METRIC[GUPXX]*METRIC_LAP_PSI4[PSIM4] - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTX]*METRIC[SHIFTX];
    g4up[1][2] = METRIC[GUPXY]*METRIC_LAP_PSI4[PSIM4] - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTX]*METRIC[SHIFTY];
    g4up[1][3] = METRIC[GUPXZ]*METRIC_LAP_PSI4[PSIM4] - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTX]*METRIC[SHIFTZ];
    g4up[2][2] = METRIC[GUPYY]*METRIC_LAP_PSI4[PSIM4] - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTY]*METRIC[SHIFTY];
    g4up[2][3] = METRIC[GUPYZ]*METRIC_LAP_PSI4[PSIM4] - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTY]*METRIC[SHIFTZ];
    g4up[3][3] = METRIC[GUPZZ]*METRIC_LAP_PSI4[PSIM4] - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTZ]*METRIC[SHIFTZ];
  }

  int count;
  // Next compute T^{\mu \nu}. This is very useful when computing fluxes and source terms in the GRMHD evolution equations.
  // (Eq. 33 in http://arxiv.org/pdf/astro-ph/0503420.pdf):
  // T^{mn} = (rho_0 h + b^2) u^m u^n + (P + 0.5 b^2) g^{mn} - b^m b^n, where m and n both run from 0 to 3.
  count=0; for(int ii=0;ii<4;ii++) for(int jj=ii;jj<4;jj++) { TUPMUNU[count] = rho0_h_plus_b2*uUP[ii]*uUP[jj] + P_plus_half_b2*g4up[ii][jj] - smallb[SMALLBT+ii]*smallb[SMALLBT+jj]; count++; }

  // Next compute T_{\mu \nu}
  // T_{mn} = (rho_0 h + b^2) u_m u_n + (P + 0.5 b^2) g_{mn} - b_m b_n, where m and n both run from 0 to 3.
  if(!already_computed_physical_metric_and_inverse) {
    CCTK_REAL LAPSE_SQUARED=SQR(METRIC_LAP_PSI4[LAPSE]);
    CCTK_REAL BETADN[4],BETAUP[4] = { 0.0, METRIC[SHIFTX],METRIC[SHIFTY],METRIC[SHIFTZ] };
    lower_4vector_output_spatial_part(METRIC_LAP_PSI4[PSI4],METRIC,BETAUP, BETADN);

    //CCTK_REAL g4dn[4][4];
    // g_{00} = - alpha^2 + gamma_{ij} beta^i beta^j = - alpha^2 beta_i beta^i
    g4dn[0][0] = -LAPSE_SQUARED + (BETAUP[1]*BETADN[1] + BETAUP[2]*BETADN[2] + BETAUP[3]*BETADN[3]);
    // g_{0i} =  gamma_{ij} beta^j = beta_i
    g4dn[0][1] = g4dn[1][0] = BETADN[1];
    g4dn[0][2] = g4dn[2][0] = BETADN[2];
    g4dn[0][3] = g4dn[3][0] = BETADN[3];
    // g_{ij} =  gamma_{ij} = Psi^4 \tilde{gamma_ij}
    g4dn[1][1] =              METRIC[GXX]*METRIC_LAP_PSI4[PSI4];
    g4dn[1][2] = g4dn[2][1] = METRIC[GXY]*METRIC_LAP_PSI4[PSI4];
    g4dn[1][3] = g4dn[3][1] = METRIC[GXZ]*METRIC_LAP_PSI4[PSI4];
    g4dn[2][2] =              METRIC[GYY]*METRIC_LAP_PSI4[PSI4];
    g4dn[2][3] = g4dn[3][2] = METRIC[GYZ]*METRIC_LAP_PSI4[PSI4];
    g4dn[3][3] =              METRIC[GZZ]*METRIC_LAP_PSI4[PSI4];
  }

  CCTK_REAL smallb_lower[NUMVARS_SMALLB];
  // b_a = b^c g_{ac}
  for(int ii=0;ii<4;ii++) { smallb_lower[SMALLBT+ii]=0; for(int jj=0;jj<4;jj++) smallb_lower[SMALLBT+ii] += smallb[SMALLBT+jj]*g4dn[ii][jj]; }

  // Compute u_0, as we've already computed u_i above.
  uDN[0]=0.0; for(int jj=0;jj<4;jj++) uDN[0] += uUP[jj]*g4dn[0][jj];

  // Compute T_{\mu \nu}
  if(update_Tmunu) {
    count=0; for(int ii=0;ii<4;ii++) for(int jj=ii;jj<4;jj++) { TDNMUNU[count] = rho0_h_plus_b2*uDN[ii]*uDN[jj] + P_plus_half_b2*g4dn[ii][jj] - smallb_lower[SMALLBT+ii]*smallb_lower[SMALLBT+jj]; count++; }
  }

  //CCTK_VInfo(CCTK_THORNSTRING,"YAY %e",TDNMUNU[0]);

  // Finally, compute conservatives:
  CONSERVS[RHOSTAR] = alpha_sqrt_gamma * U[RHOB] * uUP[0];
  CONSERVS[STILDEX] = CONSERVS[RHOSTAR]*h_enthalpy*uDN[1] + alpha_sqrt_gamma*(uUP[0]*smallb[SMALLB2]*uDN[1] - smallb[SMALLBT]*smallb_lower[SMALLBX]);
  CONSERVS[STILDEY] = CONSERVS[RHOSTAR]*h_enthalpy*uDN[2] + alpha_sqrt_gamma*(uUP[0]*smallb[SMALLB2]*uDN[2] - smallb[SMALLBT]*smallb_lower[SMALLBY]);
  CONSERVS[STILDEZ] = CONSERVS[RHOSTAR]*h_enthalpy*uDN[3] + alpha_sqrt_gamma*(uUP[0]*smallb[SMALLB2]*uDN[3] - smallb[SMALLBT]*smallb_lower[SMALLBZ]);
  CONSERVS[YET]     = CONSERVS[RHOSTAR]* U[YE];
  // tauL = alpha^2 sqrt(gamma) T^{00} - CONSERVS[RHOSTAR]
  CONSERVS[TAUENERGY] =  METRIC_LAP_PSI4[LAPSE]*alpha_sqrt_gamma*(rho0_h_plus_b2*SQR(uUP[0]) + P_plus_half_b2*(-SQR(METRIC_LAP_PSI4[LAPSEINV])) - SQR(smallb[SMALLBT])) - CONSERVS[RHOSTAR];
}
