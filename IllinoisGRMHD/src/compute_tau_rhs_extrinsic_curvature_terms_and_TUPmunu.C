static void compute_tau_rhs_extrinsic_curvature_terms_and_TUPmunu
(const int eos_key,const int keytemp,const CCTK_REAL c2pprec,const cGH *cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,CCTK_REAL *dX,CCTK_REAL **metric,gf_and_gz_struct *prims,
 CCTK_REAL **TUPmunu,eos_struct &eos,
 CCTK_REAL *gupxy,CCTK_REAL *gupxz,CCTK_REAL *gupyz,
 CCTK_REAL *kxx,CCTK_REAL *kxy,CCTK_REAL *kxz,CCTK_REAL *kyy,CCTK_REAL *kyz,CCTK_REAL *kzz,
 CCTK_REAL *tau_rhs) {

  // These loop extents must be consistent with add_fluxes_and_source_terms_to_hydro_rhss(), since we use TUPmunu there as well.
#pragma omp parallel for
  for(int k=cctk_nghostzones[2];k<cctk_lsh[2]-(cctk_nghostzones[2]-1);k++) for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-(cctk_nghostzones[1]-1);j++) for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-(cctk_nghostzones[0]-1);i++) {
        DECLARE_CCTK_PARAMETERS;

        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // First we pull in needed hydrodynamic and metric variables from memory: PART 1.
        // Reading from main memory is a SLOW operation, usually resulting in
        //   cache misses, which
        //   will waste precious runtime. Note that cache misses will often not
        //   show up when using, e.g., gprof. The slowdown due to cache misses
        //   can more than double the amount of time in this routine, so instead
        //   of reading in variables from main memory multiple times, the below
        //   forces us to read in only once, storing in local variables
        //   U{p,m}{1,2,3}, METRIC{p,m}{1,2}, U, and METRIC.

        CCTK_REAL KxxL = kxx[index];
        CCTK_REAL KxyL = kxy[index];
        CCTK_REAL KxzL = kxz[index];
        CCTK_REAL KyyL = kyy[index];
        CCTK_REAL KyzL = kyz[index];
        CCTK_REAL KzzL = kzz[index];

        //-----------------------------------------------------------------------------
        // Compute T^{\mu \nu}

        CCTK_REAL METRIC[NUMVARS_FOR_METRIC_FACEVALS]; for(int ii=0;ii<NUMVARS_FOR_METRIC_FACEVALS;ii++) METRIC[ii] = metric[ii][index];
        CCTK_REAL METRIC_LAP_PSI4[NUMVARS_METRIC_AUX]; SET_LAPSE_PSI4(METRIC_LAP_PSI4,METRIC);

        // The "vector" U represents the primitive variables: rho, P, vx, vy, vz, Bx, By, and Bz.
        CCTK_REAL U[MAXNUMVARS]; // 8 primitives in the set: {rho_b,P,vx,vy,vz,Bx,By,Bz}
        for(int ii=0;ii<MAXNUMVARS;ii++) U[ii] = prims[ii].gf[index];

        struct output_stats stats; stats.failure_checker=0;
        CCTK_REAL u0L;
        impose_speed_limit_output_u0(METRIC,U,METRIC_LAP_PSI4[PSI4],METRIC_LAP_PSI4[LAPSEINV],stats, u0L);

        /***********************************************************/
        // Compute b^{\mu} and b^2
        CCTK_REAL ONE_OVER_LAPSE = 1.0/METRIC_LAP_PSI4[LAPSE];
        CCTK_REAL ONE_OVER_LAPSE_SQRT_4PI = ONE_OVER_LAPSE*ONE_OVER_SQRT_4PI;
        CCTK_REAL u_x_over_u0_psi4,u_y_over_u0_psi4,u_z_over_u0_psi4;
        CCTK_REAL smallb[NUMVARS_SMALLB];
        compute_smallba_b2_and_u_i_over_u0_psi4(METRIC,METRIC_LAP_PSI4,U,u0L,ONE_OVER_LAPSE_SQRT_4PI,
                                                u_x_over_u0_psi4,u_y_over_u0_psi4,u_z_over_u0_psi4,smallb);

        /***********************************************************/
        // Next compute T^{\mu \nu}:
        //   First set h, the enthalpy:
        CCTK_REAL h = 0;
        if (!use_ConservativeToPrimitive){
        CCTK_REAL P_cold,eps_cold,dPcold_drho,eps_th,gamma_cold; /* <- Note that in setting h, we need to define several
                                                                        *    other variables, even though they will be unused later
                                                                        *    in this function. */
        compute_P_cold__eps_cold__dPcold_drho__eps_th__h__gamma_cold(U,eos,P_cold,eps_cold,dPcold_drho,eps_th,h,gamma_cold);
        }else{
          CCTK_REAL xrho   = U[RHOB];
          CCTK_REAL xeps   = U[EPS];
          CCTK_REAL xpress = U[PRESSURE];
          CCTK_REAL xye    = U[YE];
          CCTK_REAL xtemp  = U[TEMP];
          CCTK_INT  havetemp = 0;     // We have eps
          CCTK_INT  npoints  = 1;
          CCTK_INT  anyerr, keyerr;

          EOS_Omni_press(eos_key, keytemp, c2pprec, npoints, &xrho, &xeps, &xtemp, &xye, &xpress, &keyerr, &anyerr);

          h = 1. + xeps + xpress/xrho;
	  if(eos_key!=4){
            //for analytic EOS, recalculate h using the cold fluid quantities, this is also done in mhdflux.C and apply_tau_floor__enforce_limits_on_primitives_and_recompute_conservs.C
            CCTK_REAL P_cold,eps_cold,dPcold_drho,eps_th,gamma_cold; 
            compute_P_cold__eps_cold__dPcold_drho__eps_th__h__gamma_cold(U,eos,P_cold,eps_cold,dPcold_drho,eps_th,h,gamma_cold);
	  }
        }

        CCTK_REAL Psi6 = METRIC_LAP_PSI4[PSI2]*METRIC_LAP_PSI4[PSI4];
        CCTK_REAL Psim4 = 1.0/METRIC_LAP_PSI4[PSI4];

        CCTK_REAL rho0_h_plus_b2 = (U[RHOB]*h + smallb[SMALLB2]);
        CCTK_REAL P_plus_half_b2 = (U[PRESSURE]+0.5*smallb[SMALLB2]);

        CCTK_REAL half_alpha_sqrtgamma = 0.5*METRIC_LAP_PSI4[LAPSE]*Psi6;
        CCTK_REAL uUP[4] = { u0L, u0L*U[VX],u0L*U[VY],u0L*U[VZ] };
        // If you like, see Eq 2.119 in Numerical Relativity, by Baumgarte & Shapiro:
        CCTK_REAL ONE_OVER_LAPSE_SQUARED = SQR(ONE_OVER_LAPSE);
        CCTK_REAL g4up[4][4];

        g4up[0][0] = -ONE_OVER_LAPSE_SQUARED;
        g4up[0][1] = ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTX];
        g4up[0][2] = ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTY];
        g4up[0][3] = ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTZ];
        g4up[1][1] = METRIC[GUPXX]*Psim4 - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTX]*METRIC[SHIFTX];
        // Note that for i!=j, gupij is not stored in METRIC, since we don't need it in face value calculations.
        g4up[1][2] = gupxy[index] *Psim4 - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTX]*METRIC[SHIFTY];
        g4up[1][3] = gupxz[index] *Psim4 - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTX]*METRIC[SHIFTZ];
        g4up[2][2] = METRIC[GUPYY]*Psim4 - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTY]*METRIC[SHIFTY];
        g4up[2][3] = gupyz[index] *Psim4 - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTY]*METRIC[SHIFTZ];
        g4up[3][3] = METRIC[GUPZZ]*Psim4 - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTZ]*METRIC[SHIFTZ];

        // Next compute T^{\mu \nu}:
        // (Eq. 33 in http://arxiv.org/pdf/astro-ph/0503420.pdf):
        // T^{mn} = (rho_0 h + b^2) u^m u^n + (P + 0.5 b^2) g^{mn} - b^m b^n, where m and n both run from 0 to 3.
        CCTK_REAL TUP[4][4];
        for(int ii=0;ii<4;ii++) for(int jj=ii;jj<4;jj++) TUP[ii][jj] = rho0_h_plus_b2*uUP[ii]*uUP[jj] + P_plus_half_b2*g4up[ii][jj] - smallb[SMALLBT+ii]*smallb[SMALLBT+jj];
        //-----------------------------------------------------------------------------

        //-----------------------------------------------------------------------------
        // If we are not in the ghostzones, then compute extrinsic curvature terms for tau_rhs:
        //    Without this if() statement, _rhs variables are in general set to nonzero values in ghostzones, which messes up frozen BC's.
        //    Also, this if() statement should speed up the computation slightly.
        if(k<cctk_lsh[2]-cctk_nghostzones[2] && j<cctk_lsh[1]-cctk_nghostzones[1] && i<cctk_lsh[0]-cctk_nghostzones[0]) {
          //  \alpha \sqrt{\gamma} (T^{00} \beta^i \beta ^j + 2T^{0i} \beta^j + T^{ij})*K_{ij}
          CCTK_REAL alpha_sqrtgamma = 2.0*half_alpha_sqrtgamma;

          //  \alpha \sqrt{\gamma} (T^{00} \beta^i \beta ^j + 2T^{0i} \beta^j + T^{ij})*K_{ij}
          CCTK_REAL tau_rhs_extrinsic_curvature_terms = alpha_sqrtgamma*
            (
             TUP[0][0] * (SQR(METRIC[SHIFTX])*KxxL + SQR(METRIC[SHIFTY])*KyyL + SQR(METRIC[SHIFTZ])*KzzL +
                          2.0*(METRIC[SHIFTX]*METRIC[SHIFTY]*KxyL + METRIC[SHIFTX]*METRIC[SHIFTZ]*KxzL + METRIC[SHIFTY]*METRIC[SHIFTZ]*KyzL) ) +
             2.0*(TUP[0][1]*METRIC[SHIFTX]*KxxL + TUP[0][2]*METRIC[SHIFTY]*KyyL + TUP[0][3]*METRIC[SHIFTZ]*KzzL +
                  (TUP[0][1]*METRIC[SHIFTY] + TUP[0][2]*METRIC[SHIFTX])*KxyL +
                  (TUP[0][1]*METRIC[SHIFTZ] + TUP[0][3]*METRIC[SHIFTX])*KxzL +
                  (TUP[0][2]*METRIC[SHIFTZ] + TUP[0][3]*METRIC[SHIFTY])*KyzL ) +
             TUP[1][1]*KxxL + TUP[2][2]*KyyL + TUP[3][3]*KzzL +
             2.0*(TUP[1][2]*KxyL + TUP[1][3]*KxzL + TUP[2][3]*KyzL) );

          tau_rhs[index] = tau_rhs_extrinsic_curvature_terms;
        }

        // Set the T^{\mu \nu} gridfunction, since computing T^{\mu \nu} is expensive
        int counter=0;
        for(int ii=0;ii<4;ii++) for(int jj=ii;jj<4;jj++) { TUPmunu[counter][index] = TUP[ii][jj]; counter++; }
      }
}
