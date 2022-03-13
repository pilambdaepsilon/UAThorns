static inline CCTK_REAL fasterpow_ppm_reconstruct(CCTK_REAL inputvar,CCTK_REAL inputpow) {
  if(inputpow==2.0) return SQR(inputvar);
  return pow(inputvar,inputpow);
}

static inline void find_cp_cm(CCTK_REAL &cplus,CCTK_REAL &cminus,CCTK_REAL v02,CCTK_REAL u0,
                       CCTK_REAL vi,CCTK_REAL ONE_OVER_LAPSE_SQUARED,CCTK_REAL shifti,CCTK_REAL psim4,CCTK_REAL gupii) {
  // This computes phase speeds in the direction given by flux_dirn.
  //  Note that we replace the full dispersion relation with a simpler
  //     one, which overestimates the max. speeds by a factor of ~2.
  //     See full discussion around Eqs. 49 and 50 in
  //     http://arxiv.org/pdf/astro-ph/0503420.pdf .
  //  What follows is a complete derivation of the quadratic we solve.
  // wcm = (-k_0 u0 - k_x ux)
  // kcm^2 = K_{\mu} K^{\mu},
  // K_{\mu} K^{\mu} = (g_{\mu a} + u_{\mu} u_a) k^a * g^{\mu b} [ (g_{c b} + u_c u_b) k^c ]
  // --> g^{\mu b} (g_{c b} + u_{c} u_{b}) k^c = (\delta^{\mu}_c + u_c u^{\mu} ) k^c
  //                 = (g_{\mu a} + u_{\mu} u_a) k^a * (\delta^{\mu}_c + u_c u^{\mu} ) k^c
  //                 =[(g_{\mu a} + u_{\mu} u_a) \delta^{\mu}_c + (g_{\mu a} + u_{\mu} u_a) u_c u^{\mu} ] k^c k^a
  //                 =[(g_{c a} + u_c u_a) + (u_c u_a -  u_a u_c] k^c k^a
  //                 =(g_{c a} + u_c u_a) k^c k^a
  //                 = k_a k^a + u^c u^a k_c k_a
  // k^a = g^{\mu a} k_{\mu} = g^{0 a} k_0 + g^{x a} k_x
  // k_a k^a = k_0 g^{0 0} k_0 + k_x k_0 g^{0 x} + g^{x 0} k_0 k_x + g^{x x} k_x k_x
  //         = g^{00} (k_0)^2 + 2 g^{x0} k_0 k_x + g^{xx} (k_x)^2
  // u^c u^a k_c k_a = (u^0 k_0 + u^x k_x) (u^0 k_0 + u^x k_x) = (u^0 k_0)^2 + 2 u^x k_x u^0 k_0 + (u^x k_x)^2
  // (k_0 u0)^2  + 2 k_x ux k_0 u0 + (k_x ux)^2 = v02 [ (u^0 k_0)^2 + 2 u^x k_x u^0 k_0 + (u^x k_x)^2 + g^{00} (k_0)^2 + 2 g^{x0} k_0 k_x + g^{xx} (k_x)^2]
  // (1-v02) (u^0 k_0 + u^x k_x)^2 = v02 (g^{00} (k_0)^2 + 2 g^{x0} k_0 k_x + g^{xx} (k_x)^2)
  // (1-v02) (u^0 k_0/k_x + u^x)^2 = v02 (g^{00} (k_0/k_x)^2 + 2 g^{x0} k_0/k_x + g^{xx})
  // (1-v02) (u^0 X + u^x)^2 = v02 (g^{00} X^2 + 2 g^{x0} X + g^{xx})
  // (1-v02) (u0^2 X^2 + 2 ux u0 X + ux^2) = v02 (g^{00} X^2 + 2 g^{x0} X + g^{xx})
  // X^2 ( (1-v02) u0^2 - v02 g^{00}) + X (2 ux u0 (1-v02) - 2 v02 g^{x0}) + (1-v02) ux^2 - v02 g^{xx}
  // a = (1-v02) u0^2 - v02 g^{00} = (1-v02) u0^2 + v02/lapse^2 <-- VERIFIED
  // b = 2 ux u0 (1-v02) - 2 v02 shiftx/lapse^2 <-- VERIFIED, X->-X, because X = -w/k_1, and we are solving for -X.
  // c = (1-v02) ux^2 - v02 (gupxx*psim4 - (shiftx/lapse)^2) <-- VERIFIED
  // v02 = v_A^2 + c_s^2 (1 - v_A^2)
  CCTK_REAL u0_SQUARED=SQR(u0);

  //Find cplus, cminus:
  CCTK_REAL a = u0_SQUARED * (1.0-v02) + v02*ONE_OVER_LAPSE_SQUARED;
  CCTK_REAL b = 2.0* ( shifti*ONE_OVER_LAPSE_SQUARED * v02 - u0_SQUARED * vi * (1.0-v02) );
  CCTK_REAL c = u0_SQUARED*SQR(vi) * (1.0-v02) - v02 * ( psim4*gupii -
                                                         SQR(shifti)*ONE_OVER_LAPSE_SQUARED);
  CCTK_REAL detm = b*b - 4.0*a*c;
  //ORIGINAL LINE OF CODE:
  //if(detm < 0.0) detm = 0.0;
  //New line of code (without the if() statement) has the same effect:
  detm = sqrt(0.5*(detm + fabs(detm))); /* Based on very nice suggestion from Roland Haas */

  cplus = 0.5*(detm-b)/a;
  cminus = -0.5*(detm+b)/a;
  if (cplus < cminus) {
    CCTK_REAL cp = cminus;
    cminus = cplus;
    cplus = cp;
  }
}

static inline void compute_v02(const int eos_key,const CCTK_REAL c2pprec,CCTK_REAL dPcold_drho,CCTK_REAL gamma_therma,CCTK_REAL eps_th,CCTK_REAL h,CCTK_REAL *smallb,CCTK_REAL *U,  CCTK_REAL &v02L) {

  DECLARE_CCTK_PARAMETERS;

  if(U[RHOB]<=0) { v02L=1.0; return; }

  CCTK_REAL c_s_squared;

  if(!use_ConservativeToPrimitive){
    /* c_s = sound speed = (dP_c/drho + \Gamma(\Gamma-1) \epsilon_th)/h */
    c_s_squared = (dPcold_drho + gamma_therma*(gamma_therma-1.0)*eps_th)/(h);
  }else{
    CCTK_REAL xrho   = U[RHOB];
    CCTK_REAL xeps   = U[EPS];
    CCTK_REAL xye    = U[YE];
    CCTK_REAL xtemp  = 0.1;     //just initialization
    CCTK_INT  keytemp = 0;     // We enter with eps (reconstructed)
    CCTK_INT  npoints  = 1;
    CCTK_INT  anyerr, keyerr;

    EOS_Omni_cs2(eos_key,keytemp,c2pprec,npoints,&xrho,&xeps,&xtemp,&xye,&c_s_squared,&keyerr,&anyerr);
  }

  /* v_A = Alfven speed = b^2/(rho0 h + b^2) */
  CCTK_REAL v_A_squared = smallb[SMALLB2]/(smallb[SMALLB2] + U[RHOB]*(h));
  v02L = v_A_squared + c_s_squared*(1.0-v_A_squared);
}

static inline void compute_P_cold__eps_cold__dPcold_drho__eps_th__h__gamma_cold(CCTK_REAL *U, eos_struct &eos,
                                                                         CCTK_REAL &P_cold,CCTK_REAL &eps_cold,CCTK_REAL &dPcold_drho,CCTK_REAL &eps_th,CCTK_REAL &h,
                                                                         CCTK_REAL &gamma_cold) {
  // This code handles equations of state of the form defined
  // in Eqs 13-16 in http://arxiv.org/pdf/0802.0200.pdf

  if(U[RHOB]==0) {
    P_cold = 0.0;
    eps_cold = 0.0;
    dPcold_drho = 0.0;
    eps_th = 0.0;
    h = 0.0;
    gamma_cold = eos.gamma_tab[0];
    return;
  }

  CCTK_REAL U_RHOB_inv = 1.0/U[RHOB];

  if(eos.neos==1) {
    // Eq. 14 of http://arxiv.org/pdf/0802.0200.pdf :
    // P_{cold} = K_i rho_i^{\Gamma_i}
    P_cold = eos.k_tab[0]*fasterpow_ppm_reconstruct(U[RHOB],eos.gamma_tab[0]);
    // Eq. 16 of http://arxiv.org/pdf/0802.0200.pdf :
    // \epsilon_{cold} = \int ( P_{cold}(rho) / rho^2 ) drho
    //                 = \int ( K_0 \rho^{\Gamma_0 - 2} ) drho
    //                 = ( K_0 \rho^{\Gamma_0 - 1} ) / (\Gamma_0 - 1)
    //                 = ( P_{cold} / rho          ) / (\Gamma_0 - 1)
    eps_cold = P_cold*U_RHOB_inv/(eos.gamma_tab[0]-1.0);
    // dPcold/drho = K_i \Gamma_i rho_i^{\Gamma_i-1} = \Gamma_i P_{cold} / rho
    dPcold_drho = eos.gamma_tab[0]*P_cold*U_RHOB_inv;
    // Eq. 15 of http://arxiv.org/pdf/0802.0200.pdf :
    // P_{th} = (\Gamma_{th} - 1) \rho_0 \epsilon_{th},
    // Eq. 13 of http://arxiv.org/pdf/0802.0200.pdf :
    // P_{th} = P - P_{cold}
    // -> P - P_{cold} = (\Gamma_{th} - 1) \rho_0 \epsilon_{th}
    // -> \epsilon_{th} = ( P - P_{cold} ) / [ (\Gamma_{th} - 1) \rho_0 ]
    eps_th = (U[PRESSURE] - P_cold)/(eos.gamma_th-1.0)*U_RHOB_inv;
    // Just below Eq. 16 in http://arxiv.org/pdf/astro-ph/0503420.pdf :
    // h = 1 + \epsilon + P/rho
    h = 1.0 + eps_cold + eps_th + U[PRESSURE]*U_RHOB_inv;
    gamma_cold = eos.gamma_tab[0];
    return;
  }

  // See comments above for the eos.neos==1 case for relevant
  // equations & references; the extension to arbitrary "nn"
  // is straightforward.
  for(int nn=1;nn<eos.neos;nn++) {
    if (U[RHOB] <= eos.rho_tab[nn]  &&  U[RHOB] > eos.rho_tab[nn-1]) {
      P_cold = eos.k_tab[nn]*fasterpow_ppm_reconstruct(U[RHOB],eos.gamma_tab[nn]);
      eps_cold = eos.eps_tab[nn-1] + (P_cold*U_RHOB_inv - eos.P_tab[nn-1]/eos.rho_tab[nn-1])/(eos.gamma_tab[nn]-1.0);
      dPcold_drho = eos.gamma_tab[nn]*P_cold*U_RHOB_inv;
      eps_th = (U[PRESSURE] - P_cold)/(eos.gamma_th-1.0)*U_RHOB_inv;
      h = 1.0 + eps_cold + eps_th + U[PRESSURE]*U_RHOB_inv;
      gamma_cold = eos.gamma_tab[nn];
    }
  }
  if (U[RHOB] > eos.rho_tab[eos.neos-1]) {
    P_cold = eos.k_tab[eos.neos]*fasterpow_ppm_reconstruct(U[RHOB],eos.gamma_tab[eos.neos]);
    eps_cold = eos.eps_tab[eos.neos-1] + (P_cold*U_RHOB_inv - eos.P_tab[eos.neos-1]/eos.rho_tab[eos.neos-1])/(eos.gamma_tab[eos.neos]-1.0);
    dPcold_drho = eos.gamma_tab[eos.neos]*P_cold*U_RHOB_inv;
    eps_th = (U[PRESSURE] - P_cold)/(eos.gamma_th-1.0)*U_RHOB_inv;
    h = 1.0 + eps_cold + eps_th + U[PRESSURE]*U_RHOB_inv;
    gamma_cold = eos.gamma_tab[eos.neos];
  }
}

// b_x = g_{\mu x} b^{\mu}
//     = g_{t x} b^t + g_{i x} b^i
//     = b^t gamma_{xj} beta^j + gamma_{ix} b^i
//     = gamma_{xj} (b^j + beta^j b^t)
static inline void lower_4vector_output_spatial_part(CCTK_REAL psi4,CCTK_REAL *METRIC,CCTK_REAL *smallb, CCTK_REAL *smallb_lower) {
  smallb_lower[SMALLBX] = psi4*( METRIC[GXX]*(smallb[SMALLBX]+smallb[SMALLBT]*METRIC[SHIFTX]) + METRIC[GXY]*(smallb[SMALLBY]+smallb[SMALLBT]*METRIC[SHIFTY]) +
                                 METRIC[GXZ]*(smallb[SMALLBZ]+smallb[SMALLBT]*METRIC[SHIFTZ]) );
  smallb_lower[SMALLBY] = psi4*( METRIC[GXY]*(smallb[SMALLBX]+smallb[SMALLBT]*METRIC[SHIFTX]) + METRIC[GYY]*(smallb[SMALLBY]+smallb[SMALLBT]*METRIC[SHIFTY]) +
                                 METRIC[GYZ]*(smallb[SMALLBZ]+smallb[SMALLBT]*METRIC[SHIFTZ]) );
  smallb_lower[SMALLBZ] = psi4*( METRIC[GXZ]*(smallb[SMALLBX]+smallb[SMALLBT]*METRIC[SHIFTX]) + METRIC[GYZ]*(smallb[SMALLBY]+smallb[SMALLBT]*METRIC[SHIFTY]) +
                                 METRIC[GZZ]*(smallb[SMALLBZ]+smallb[SMALLBT]*METRIC[SHIFTZ]) );
}

static inline void impose_speed_limit_output_u0(CCTK_REAL *METRIC,CCTK_REAL *U,CCTK_REAL psi4,CCTK_REAL ONE_OVER_LAPSE,output_stats &stats, CCTK_REAL &u0_out) {
  DECLARE_CCTK_PARAMETERS;

  // Derivation of first equation:
  // \gamma_{ij} (v^i + \beta^i)(v^j + \beta^j)/(\alpha)^2
  //   = \gamma_{ij} 1/(u^0)^2 ( \gamma^{ik} u_k \gamma^{jl} u_l /(\alpha)^2 <- Using Eq. 53 of arXiv:astro-ph/0503420
  //   = 1/(u^0 \alpha)^2 u_j u_l \gamma^{jl}  <- Since \gamma_{ij} \gamma^{ik} = \delta^k_j
  //   = 1/(u^0 \alpha)^2 ( (u^0 \alpha)^2 - 1 ) <- Using Eq. 56 of arXiv:astro-ph/0503420
  //   = 1 - 1/(u^0 \alpha)^2 <= 1
  CCTK_REAL one_minus_one_over_alpha_u0_squared = psi4*(METRIC[GXX]* SQR(U[VX] + METRIC[SHIFTX]) +
                                                        2.0*METRIC[GXY]*(U[VX] + METRIC[SHIFTX])*(U[VY] + METRIC[SHIFTY]) +
                                                        2.0*METRIC[GXZ]*(U[VX] + METRIC[SHIFTX])*(U[VZ] + METRIC[SHIFTZ]) +
                                                        METRIC[GYY]* SQR(U[VY] + METRIC[SHIFTY]) +
                                                        2.0*METRIC[GYZ]*(U[VY] + METRIC[SHIFTY])*(U[VZ] + METRIC[SHIFTZ]) +
                                                        METRIC[GZZ]* SQR(U[VZ] + METRIC[SHIFTZ]) )*SQR(ONE_OVER_LAPSE);

  /*** Limit velocity to GAMMA_SPEED_LIMIT ***/
  const CCTK_REAL ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED = 1.0-1.0/SQR(GAMMA_SPEED_LIMIT);
  if(one_minus_one_over_alpha_u0_squared > ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED) {
    CCTK_REAL correction_fac = sqrt(ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED/one_minus_one_over_alpha_u0_squared);
    U[VX] = (U[VX] + METRIC[SHIFTX])*correction_fac-METRIC[SHIFTX];
    U[VY] = (U[VY] + METRIC[SHIFTY])*correction_fac-METRIC[SHIFTY];
    U[VZ] = (U[VZ] + METRIC[SHIFTZ])*correction_fac-METRIC[SHIFTZ];
    one_minus_one_over_alpha_u0_squared=ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED;
    stats.failure_checker+=1000;
  }

  // A = 1.0-one_minus_one_over_alpha_u0_squared = 1-(1-1/(al u0)^2) = 1/(al u0)^2
  // 1/sqrt(A) = al u0
  //CCTK_REAL alpha_u0_minus_one = 1.0/sqrt(1.0-one_minus_one_over_alpha_u0_squared)-1.0;
  //u0_out          = (alpha_u0_minus_one + 1.0)*ONE_OVER_LAPSE;
  CCTK_REAL alpha_u0 = 1.0/sqrt(1.0-one_minus_one_over_alpha_u0_squared);
  if(std::isnan(alpha_u0*ONE_OVER_LAPSE)) printf("BAD FOUND NAN U0 CALC: %.15e %.15e %.15e | %.15e %.15e\n",alpha_u0,ONE_OVER_LAPSE,one_minus_one_over_alpha_u0_squared,psi4, U[VX]);
  u0_out = alpha_u0*ONE_OVER_LAPSE;
}

// The two lines of code below are written to reduce roundoff error and were in the above function. I don't think they reduce error.
// one_over_alpha_u0  = sqrt(1.0-one_minus_one_over_alpha_u0_squared);
/* Proof of following line: */
/*   [        1-1/(alphau0)^2         ] / [ 1/(alphau0) (1 + 1/(alphau0)) ] */
/* = [ (alphau0)^2 - 1)/((alphau0)^2) ] / [ 1/(alphau0) + 1/(alphau0)^2   ] */
/* = [ (alphau0)^2 - 1)/((alphau0)^2) ] / [  (alphau0 + 1)/(alphau0)^2    ] */
/* = [        (alphau0)^2 - 1)        ] / [        (alphau0 + 1)          ] */
/*   [   (alphau0 + 1) (alphau0 - 1)  ] / [        (alphau0 + 1)          ] */
/* =                              alphau0 - 1                               */
//alpha_u0_minus_one = one_minus_one_over_alpha_u0_squared/one_over_alpha_u0/(1.0+one_over_alpha_u0);
//u0_out = (alpha_u0_minus_one+1.0)*ONE_OVER_LAPSE;

static inline void enforce_pressure_floor_ceiling(output_stats &stats,CCTK_REAL kpoly,CCTK_REAL P_cold,CCTK_REAL Psi6,const CCTK_REAL Psi6threshold,CCTK_REAL rho_b,const CCTK_REAL rhobatm,  CCTK_REAL &P) {
  CCTK_REAL P_min=0.9*P_cold;
  if(P<P_min) {
    stats.failure_checker+=10;
    P=P_min;
  }
  //MAX(P,P_min);
  //if(P < P_min) P=1.0*P_cold;

  /* OLD: Discarded because lower limit is unphysical.
     if(P <= 0.5*kpoly*P_cold) {
     P=0.5*kpoly*P_cold;
     }
  */

  //CCTK_REAL P_max = 10.0*P_cold;
  CCTK_REAL P_max = 100.0*P_cold;
  if(Psi6 > Psi6threshold) P_max = 1e5*P_cold; // <-- better than 10.

  if((rho_b < 100.0*rhobatm || Psi6 > Psi6threshold) && P>P_max) {
    P=P_max;
    stats.failure_checker+=100;
  }

  /*
    CCTK_REAL rho_horiz_cap = 1000.0*rhobatm;

    //New density damping mechanism inside the horizon
    if(Psi6 > Psi6threshold && rho_b>rho_horiz_cap) {
    CCTK_REAL six_phi=log(Psi6);
    CCTK_REAL six_phithreshold=log(Psi6threshold);
    CCTK_REAL Psi6max_approx=350000;
    rho_b = rho_horiz_cap+(rho_b-rho_horiz_cap)*exp(-200.0*SQR((six_phi-six_phithreshold)/log(Psi6max_approx)));
    }
  */
}

static inline void compute_smallba_b2_and_u_i_over_u0_psi4(CCTK_REAL *METRIC,CCTK_REAL *METRIC_LAP_PSI4,CCTK_REAL *U,CCTK_REAL u0L,CCTK_REAL ONE_OVER_LAPSE_SQRT_4PI,
                                                    CCTK_REAL &u_x_over_u0_psi4,CCTK_REAL &u_y_over_u0_psi4,CCTK_REAL &u_z_over_u0_psi4,CCTK_REAL *smallb) {

  // NOW COMPUTE b^{\mu} and b^2 = b^{\mu} b^{\nu} g_{\mu \nu}
  CCTK_REAL ONE_OVER_U0 = 1.0/u0L;
  CCTK_REAL shiftx_plus_vx = (METRIC[SHIFTX]+U[VX]);
  CCTK_REAL shifty_plus_vy = (METRIC[SHIFTY]+U[VY]);
  CCTK_REAL shiftz_plus_vz = (METRIC[SHIFTZ]+U[VZ]);

  // Eq. 56 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
  //  u_i = gamma_{ij} u^0 (v^j + beta^j), gamma_{ij} is the physical metric, and gamma_{ij} = Psi4 * METRIC[Gij], since METRIC[Gij] is the conformal metric.
  u_x_over_u0_psi4 =  METRIC[GXX]*shiftx_plus_vx + METRIC[GXY]*shifty_plus_vy + METRIC[GXZ]*shiftz_plus_vz;
  u_y_over_u0_psi4 =  METRIC[GXY]*shiftx_plus_vx + METRIC[GYY]*shifty_plus_vy + METRIC[GYZ]*shiftz_plus_vz;
  u_z_over_u0_psi4 =  METRIC[GXZ]*shiftx_plus_vx + METRIC[GYZ]*shifty_plus_vy + METRIC[GZZ]*shiftz_plus_vz;

  // Eqs. 23 and 31 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
  //   Compute alpha sqrt(4 pi) b^t = u_i B^i
  CCTK_REAL alpha_sqrt_4pi_bt = ( u_x_over_u0_psi4*U[BX_CENTER] + u_y_over_u0_psi4*U[BY_CENTER] + u_z_over_u0_psi4*U[BZ_CENTER] ) * METRIC_LAP_PSI4[PSI4]*u0L;
  // Eq. 24 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
  // b^i = B^i_u / sqrt(4 pi)
  // b^i = ( B^i/alpha + B^0_u u^i ) / ( u^0 sqrt(4 pi) )
  // b^i = ( B^i/alpha +  sqrt(4 pi) b^t u^i ) / ( u^0 sqrt(4 pi) )
  // b^i = ( B^i +  alpha sqrt(4 pi) b^t u^i ) / ( alpha u^0 sqrt(4 pi) )
  // b^i = ( B^i/u^0 +  alpha sqrt(4 pi) b^t u^i/u^0 ) / ( alpha sqrt(4 pi) )
  // b^i = ( B^i/u^0 +  alpha sqrt(4 pi) b^t v^i ) / ( alpha sqrt(4 pi) )
  smallb[SMALLBX] = (U[BX_CENTER]*ONE_OVER_U0 + U[VX]*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
  smallb[SMALLBY] = (U[BY_CENTER]*ONE_OVER_U0 + U[VY]*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
  smallb[SMALLBZ] = (U[BZ_CENTER]*ONE_OVER_U0 + U[VZ]*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
  // Eq. 23 in http://arxiv.org/pdf/astro-ph/0503420.pdf, with alpha sqrt (4 pi) b^2 = u_i B^i already computed above
  smallb[SMALLBT] = alpha_sqrt_4pi_bt * ONE_OVER_LAPSE_SQRT_4PI;

  // b^2 = g_{\mu \nu} b^{\mu} b^{\nu}
  //     = gtt bt^2 + gxx bx^2 + gyy by^2 + gzz bz^2 + 2 (gtx bt bx + gty bt by + gtz bt bz + gxy bx by + gxz bx bz + gyz by bz)
  //     = (-al^2 + gamma_{ij} betai betaj) bt^2 + b^i b^j gamma_{ij} + 2 g_{t i} b^t b^i
  //     = - (alpha b^t)^2 + (b^t)^2 gamma_{ij} beta^i beta^j + b^i b^j gamma_{ij} + 2 b^t g_{t i} b^i
  //     = - (alpha b^t)^2 + (b^t)^2 gamma_{ij} beta^i beta^j + b^i b^j gamma_{ij} + 2 b^t (gamma_{ij} beta^j) b^i
  //     = - (alpha b^t)^2 + gamma_{ij} ((b^t)^2 beta^i beta^j + b^i b^j + 2 b^t beta^j b^i)
  //     = - (alpha b^t)^2 + gamma_{ij} ((b^t)^2 beta^i beta^j + 2 b^t beta^j b^i + b^i b^j)
  //     = - (alpha b^t)^2 + gamma_{ij} (b^i + b^t beta^i) (b^j + b^t beta^j)
  CCTK_REAL bx_plus_shiftx_bt = smallb[SMALLBX]+METRIC[SHIFTX]*smallb[SMALLBT];
  CCTK_REAL by_plus_shifty_bt = smallb[SMALLBY]+METRIC[SHIFTY]*smallb[SMALLBT];
  CCTK_REAL bz_plus_shiftz_bt = smallb[SMALLBZ]+METRIC[SHIFTZ]*smallb[SMALLBT];
  smallb[SMALLB2] = -SQR(METRIC_LAP_PSI4[LAPSE]*smallb[SMALLBT]) +
    (  METRIC[GXX]*SQR(bx_plus_shiftx_bt) + METRIC[GYY]*SQR(by_plus_shifty_bt) + METRIC[GZZ]*SQR(bz_plus_shiftz_bt) +
       2.0*( METRIC[GXY]*(bx_plus_shiftx_bt)*(by_plus_shifty_bt) +
             METRIC[GXZ]*(bx_plus_shiftx_bt)*(bz_plus_shiftz_bt) +
             METRIC[GYZ]*(by_plus_shifty_bt)*(bz_plus_shiftz_bt) )  ) * METRIC_LAP_PSI4[PSI4]; // mult by psi4 because METRIC[GIJ] is the conformal metric.
  /***********************************************************/
}
