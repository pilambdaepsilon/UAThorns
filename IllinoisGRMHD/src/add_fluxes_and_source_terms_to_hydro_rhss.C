// Side note: the following values could be used for cell averaged gfs:
//     am2=-1.0/12.0, am1=7.0/12.0, a0=7.0/12.0, a1=-1.0/12.0
// However, since the metric gfs store the grid point values instead of the cell average,
//     the following coefficients should be used:
//     am2 = -1/16, am1 = 9/16, a0 = 9/16, a1 = -1/16
// This will yield the third-order-accurate face values at m-1/2,
//      using values specified at {m-2,m-1,m,m+1}
#define AM2 -0.0625
#define AM1  0.5625
#define A0   0.5625
#define A1  -0.0625
#define COMPUTE_FCVAL(METRICm2,METRICm1,METRIC,METRICp1) (AM2*(METRICm2) + AM1*(METRICm1) + A0*(METRIC) + A1*(METRICp1))

static inline void mhdflux(const int eos_key,const CCTK_REAL c2pprec,int i,int j,int k,const int flux_dirn,CCTK_REAL *Ul,CCTK_REAL *Ur,  CCTK_REAL *FACEVAL,CCTK_REAL *FACEVAL_LAPSE_PSI4,eos_struct &eos,
                           CCTK_REAL &cmax,CCTK_REAL &cmin,
                           CCTK_REAL &rho_star_flux,CCTK_REAL &tau_flux,CCTK_REAL &st_x_flux,CCTK_REAL &st_y_flux,CCTK_REAL &st_z_flux, CCTK_REAL &Yet_flux);

#define COMPUTE_FOURMETRIC(g4tt,g4tx,g4ty,g4tz,g4xx,g4xy,g4xz,g4yy,g4yz,g4zz,METRIC,METRIC_AUX)  ( { \
      /* g_{0i} = beta_i */                                             \
      g4tx = METRIC_AUX[PSI4]*(METRIC[GXX]*METRIC[SHIFTX] + METRIC[GXY]*METRIC[SHIFTY] + METRIC[GXZ]*METRIC[SHIFTZ]); \
      g4ty = METRIC_AUX[PSI4]*(METRIC[GXY]*METRIC[SHIFTX] + METRIC[GYY]*METRIC[SHIFTY] + METRIC[GYZ]*METRIC[SHIFTZ]); \
      g4tz = METRIC_AUX[PSI4]*(METRIC[GXZ]*METRIC[SHIFTX] + METRIC[GYZ]*METRIC[SHIFTY] + METRIC[GZZ]*METRIC[SHIFTZ]); \
      /* g_{00} = -alpha^2 + beta^i beta^j gamma_{ij} = -alpha^2 + beta^i beta_i = -alpha^2 + beta^i g_{0i} */ \
      g4tt = -SQR(METRIC_AUX[LAPSE]) + g4tx*METRIC[SHIFTX] + g4ty*METRIC[SHIFTY] + g4tz*METRIC[SHIFTZ]; \
      g4xx = METRIC_AUX[PSI4]*METRIC[GXX];                              \
      g4xy = METRIC_AUX[PSI4]*METRIC[GXY];                              \
      g4xz = METRIC_AUX[PSI4]*METRIC[GXZ];                              \
      g4yy = METRIC_AUX[PSI4]*METRIC[GYY];                              \
      g4yz = METRIC_AUX[PSI4]*METRIC[GYZ];                              \
      g4zz = METRIC_AUX[PSI4]*METRIC[GZZ];                              \
    } )


static void add_fluxes_and_source_terms_to_hydro_rhss(const int eos_key,const int keytemp,const CCTK_REAL c2pprec,const int flux_dirn,const cGH *cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,CCTK_REAL *dX,
                                                      CCTK_REAL **metric,gf_and_gz_struct *in_prims,CCTK_REAL **TUPmunu,
                                                      int numvars_reconstructed,const int *which_prims_reconstructed,gf_and_gz_struct *out_prims_r,gf_and_gz_struct *out_prims_l,eos_struct &eos,
                                                      CCTK_REAL *cmax,CCTK_REAL *cmin,
                                                      CCTK_REAL *rho_star_flux,CCTK_REAL *tau_flux,CCTK_REAL *st_x_flux,CCTK_REAL *st_y_flux,CCTK_REAL *st_z_flux, CCTK_REAL *Yet_flux,
                                                      CCTK_REAL *rho_star_rhs,CCTK_REAL *tau_rhs,CCTK_REAL *st_x_rhs,CCTK_REAL *st_y_rhs,CCTK_REAL *st_z_rhs, CCTK_REAL *Yet_rhs) {

  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL dxi[4] = { 1e100,1.0/dX[0],1.0/dX[1],1.0/dX[2] };

  // Notice in the loop below that we go from 3 to cctk_lsh-2 for i, j, AND k, even though
  //   we are only computing the flux in one direction at a time. This is because in the end,
  //   we only need the rhs's from 3 to cctk_lsh-3 for i, j, and k.
#pragma omp parallel for
  for(int k=cctk_nghostzones[2];k<cctk_lsh[2]-(cctk_nghostzones[2]-1);k++) for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-(cctk_nghostzones[1]-1);j++) for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-(cctk_nghostzones[0]-1);i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	// Set metric and associated variables
	CCTK_REAL METRIC[NUMVARS_FOR_METRIC_FACEVALS]; for(int ii=0;ii<NUMVARS_FOR_METRIC_FACEVALS;ii++) METRIC[ii] = metric[ii][index];
	CCTK_REAL METRIC_LAP_PSI4[NUMVARS_METRIC_AUX]; SET_LAPSE_PSI4(METRIC_LAP_PSI4,METRIC);

	CCTK_REAL Ur[MAXNUMVARS];
//        for(int ii=0;ii<numvars_reconstructed;ii++) Ur[which_prims_reconstructed[ii]] = out_prims_r[which_prims_reconstructed[ii]].gf[index];
	for(int ii=0;ii<numvars_reconstructed;ii++){Ur[ii] = out_prims_r[ii].gf[index];}Ur[YE]=out_prims_r[YE].gf[index];Ur[EPS]=out_prims_r[EPS].gf[index];Ur[TEMP]=out_prims_r[TEMP].gf[index];

	CCTK_REAL Ul[MAXNUMVARS];
//        for(int ii=0;ii<numvars_reconstructed;ii++) Ul[which_prims_reconstructed[ii]] = out_prims_l[which_prims_reconstructed[ii]].gf[index];
	for(int ii=0;ii<numvars_reconstructed;ii++){Ul[ii] = out_prims_l[ii].gf[index];}Ul[YE]=out_prims_l[YE].gf[index];Ul[EPS]=out_prims_l[EPS].gf[index];Ul[TEMP]=out_prims_l[TEMP].gf[index];
//	CCTK_REAL Ur[MAXNUMVARS]; for(int ii=0;ii<numvars_reconstructed;ii++) Ur[which_prims_reconstructed[ii]] = out_prims_r[which_prims_reconstructed[ii]].gf[index];
//	CCTK_REAL Ul[MAXNUMVARS]; for(int ii=0;ii<numvars_reconstructed;ii++) Ul[which_prims_reconstructed[ii]] = out_prims_l[which_prims_reconstructed[ii]].gf[index];

	// Read the T^{\mu \nu} gridfunction from memory, since computing T^{\mu \nu} is expensive
	CCTK_REAL TUP[4][4]; int counter=0;
	for(int ii=0;ii<4;ii++) for(int jj=ii;jj<4;jj++) { TUP[ii][jj] = TUPmunu[counter][index]; counter++; }

	// Next set metric on the faces, applying a 3rd-order lopsided stencil.
	int indexm2 = CCTK_GFINDEX3D(cctkGH,i-2*kronecker_delta[flux_dirn][0],j-2*kronecker_delta[flux_dirn][1],k-2*kronecker_delta[flux_dirn][2]);
	int indexm1 = CCTK_GFINDEX3D(cctkGH,i-  kronecker_delta[flux_dirn][0],j-  kronecker_delta[flux_dirn][1],k-  kronecker_delta[flux_dirn][2]);
	int indexp1 = CCTK_GFINDEX3D(cctkGH,i+  kronecker_delta[flux_dirn][0],j+  kronecker_delta[flux_dirn][1],k+  kronecker_delta[flux_dirn][2]);
	int indexp2 = CCTK_GFINDEX3D(cctkGH,i+2*kronecker_delta[flux_dirn][0],j+2*kronecker_delta[flux_dirn][1],k+2*kronecker_delta[flux_dirn][2]);
	// The "vector" METRIC stores needed metric-related quantities.
	CCTK_REAL METRICm2[NUMVARS_FOR_METRIC_FACEVALS]; for(int ii=0;ii<NUMVARS_FOR_METRIC_FACEVALS;ii++) METRICm2[ii] = metric[ii][indexm2];
	CCTK_REAL METRICm1[NUMVARS_FOR_METRIC_FACEVALS]; for(int ii=0;ii<NUMVARS_FOR_METRIC_FACEVALS;ii++) METRICm1[ii] = metric[ii][indexm1];
	CCTK_REAL METRICp1[NUMVARS_FOR_METRIC_FACEVALS]; for(int ii=0;ii<NUMVARS_FOR_METRIC_FACEVALS;ii++) METRICp1[ii] = metric[ii][indexp1];
	CCTK_REAL METRICp2[NUMVARS_FOR_METRIC_FACEVALS]; for(int ii=0;ii<NUMVARS_FOR_METRIC_FACEVALS;ii++) METRICp2[ii] = metric[ii][indexp2];

	// Next compute the metric values at the {i,j,k} +/- 1/2 faces (i.e., the "face values" of the metric)
	CCTK_REAL FACEVAL[NUMVARS_FOR_METRIC_FACEVALS],FACEVALp1[NUMVARS_FOR_METRIC_FACEVALS];
	for(int w=0;w<NUMVARS_FOR_METRIC_FACEVALS;w++) FACEVAL[w]   = COMPUTE_FCVAL(METRICm2[w],METRICm1[w],METRIC[w],METRICp1[w]);
	for(int w=0;w<NUMVARS_FOR_METRIC_FACEVALS;w++) FACEVALp1[w] = COMPUTE_FCVAL(METRICm1[w],METRIC[w],METRICp1[w],METRICp2[w]);
	// Then compute the lapse and Psi4 = exp(4*phi)
	CCTK_REAL FACEVAL_LAPSE_PSI4[NUMVARS_METRIC_AUX],FACEVAL_LAPSE_PSI4p1[NUMVARS_METRIC_AUX];
	SET_LAPSE_PSI4(FACEVAL_LAPSE_PSI4,  FACEVAL);
	SET_LAPSE_PSI4(FACEVAL_LAPSE_PSI4p1,FACEVALp1);

	//-----------------------------------------------------------------------------
	// Next compute fluxes for \tilde{S}_i, tau, and rho_*
	mhdflux(eos_key,c2pprec,i,j,k,flux_dirn,Ul  ,Ur  ,FACEVAL  ,FACEVAL_LAPSE_PSI4  ,eos, cmax[index],cmin[index],
		rho_star_flux[index],tau_flux[index],st_x_flux[index],st_y_flux[index],st_z_flux[index], Yet_flux[index]);


	//-----------------------------------------------------------------------------
	// If we are not in the ghostzones, then add third-order accurate curvature terms to \tilde{S}_i RHS's
	//    Without this if() statement, _rhs variables are in general set to nonzero values in ghostzones, which messes up frozen BC's.
	//    Also, this if() statement should speed up the computation slightly.
	if(k<cctk_lsh[2]-cctk_nghostzones[2] && j<cctk_lsh[1]-cctk_nghostzones[1] && i<cctk_lsh[0]-cctk_nghostzones[0]) {

	  CCTK_REAL Psi6 = METRIC_LAP_PSI4[PSI2]*METRIC_LAP_PSI4[PSI4];
	  CCTK_REAL half_alpha_sqrtgamma = 0.5*METRIC_LAP_PSI4[LAPSE]*Psi6;

	  // First compute four metric.
	  CCTK_REAL g4tt_f,g4tx_f,g4ty_f,g4tz_f,g4xx_f,g4xy_f,g4xz_f,g4yy_f,g4yz_f,g4zz_f;
	  COMPUTE_FOURMETRIC(g4tt_f,g4tx_f,g4ty_f,g4tz_f,g4xx_f,g4xy_f,g4xz_f,g4yy_f,g4yz_f,g4zz_f,FACEVAL,FACEVAL_LAPSE_PSI4);

	  CCTK_REAL g4tt_fp1,g4tx_fp1,g4ty_fp1,g4tz_fp1,g4xx_fp1,g4xy_fp1,g4xz_fp1,g4yy_fp1,g4yz_fp1,g4zz_fp1;
	  COMPUTE_FOURMETRIC(g4tt_fp1,g4tx_fp1,g4ty_fp1,g4tz_fp1,g4xx_fp1,g4xy_fp1,g4xz_fp1,g4yy_fp1,g4yz_fp1,g4zz_fp1,FACEVALp1,FACEVAL_LAPSE_PSI4p1);

	  // Compute \partial_i g_{\mu \nu} at m+1/2
	  CCTK_REAL partial_i_gmunu[4][4];
	  partial_i_gmunu[0][0] = (g4tt_fp1 - g4tt_f)*dxi[flux_dirn];
	  partial_i_gmunu[0][1] = (g4tx_fp1 - g4tx_f)*dxi[flux_dirn];
	  partial_i_gmunu[0][2] = (g4ty_fp1 - g4ty_f)*dxi[flux_dirn];
	  partial_i_gmunu[0][3] = (g4tz_fp1 - g4tz_f)*dxi[flux_dirn];
	  partial_i_gmunu[1][1] = (g4xx_fp1 - g4xx_f)*dxi[flux_dirn];
	  partial_i_gmunu[1][2] = (g4xy_fp1 - g4xy_f)*dxi[flux_dirn];
	  partial_i_gmunu[1][3] = (g4xz_fp1 - g4xz_f)*dxi[flux_dirn];
	  partial_i_gmunu[2][2] = (g4yy_fp1 - g4yy_f)*dxi[flux_dirn];
	  partial_i_gmunu[2][3] = (g4yz_fp1 - g4yz_f)*dxi[flux_dirn];
	  partial_i_gmunu[3][3] = (g4zz_fp1 - g4zz_f)*dxi[flux_dirn];

	  // Needed for tau_rhs computation:
	  CCTK_REAL lapse_deriv[4] = { 0,0,0,0 };
	  lapse_deriv[flux_dirn] = (FACEVALp1[LAPM1] - FACEVAL[LAPM1])*dxi[flux_dirn];

	  // Needed for st_i_rhs computation:
	  CCTK_REAL st_i_curvature_terms[4] = { 0,0,0,0 };
	  // add \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu} \partial_i g_{\mu \nu} . Note that i is given by the flux direction.
	  //   (Source term of Eq 43 in http://arxiv.org/pdf/astro-ph/0503420.pdf)
	  st_i_curvature_terms[flux_dirn] = half_alpha_sqrtgamma * ( TUP[0][0]*partial_i_gmunu[0][0] +
								     TUP[1][1]*partial_i_gmunu[1][1] +
								     TUP[2][2]*partial_i_gmunu[2][2] +
								     TUP[3][3]*partial_i_gmunu[3][3] +
								     2.0*(TUP[0][1]*partial_i_gmunu[0][1] +
									  TUP[0][2]*partial_i_gmunu[0][2] +
									  TUP[0][3]*partial_i_gmunu[0][3] +
									  TUP[1][2]*partial_i_gmunu[1][2] +
									  TUP[1][3]*partial_i_gmunu[1][3] +
									  TUP[2][3]*partial_i_gmunu[2][3]) );

	  // add - ( T^{00} \beta^i + T^{0i} ) \partial_i \alpha.
	  //   (Last part of Eq. 39 source term in http://arxiv.org/pdf/astro-ph/0503420.pdf)
	  CCTK_REAL alpha_sqrtgamma = 2.0*half_alpha_sqrtgamma;
	  tau_rhs[index]  += alpha_sqrtgamma*(-(TUP[0][0]*METRIC[SHIFTX+(flux_dirn-1)] + TUP[0][flux_dirn])*lapse_deriv[flux_dirn]);

	  // Eq 43 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
	  // \partial_t \tilde{S}_i = - \partial_i (\alpha \sqrt{\gamma} T^j_i) + \frac{1}{2}\alpha \sqrt{\gamma} T^{\mu \nu}g_{\mu \nu,i}
	  // Notice that st_i_curvature_terms[N]=0 for N!=flux_dirn.
	  st_x_rhs[index] += st_i_curvature_terms[1];
	  st_y_rhs[index] += st_i_curvature_terms[2];
	  st_z_rhs[index] += st_i_curvature_terms[3];
	}

      }

  // Notice in the loop below that we go from 3 to cctk_lsh-3 for i, j, AND k, even though
  //   we are only computing the flux in one direction. This is because in the end,
  //   we only need the rhs's from 3 to cctk_lsh-3 for i, j, and k.
#pragma omp parallel for
  for(int k=cctk_nghostzones[2];k<cctk_lsh[2]-cctk_nghostzones[2];k++) for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-cctk_nghostzones[1];j++) for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-cctk_nghostzones[0];i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	int indexp1 = CCTK_GFINDEX3D(cctkGH,i+kronecker_delta[flux_dirn][0],j+kronecker_delta[flux_dirn][1],k+kronecker_delta[flux_dirn][2]);

	rho_star_rhs[index] += (rho_star_flux[index] - rho_star_flux[indexp1]) * dxi[flux_dirn];
	tau_rhs[index]      += (tau_flux[index]      - tau_flux[indexp1]     ) * dxi[flux_dirn];
	st_x_rhs[index]     += (st_x_flux[index]     - st_x_flux[indexp1]    ) * dxi[flux_dirn];
	st_y_rhs[index]     += (st_y_flux[index]     - st_y_flux[indexp1]    ) * dxi[flux_dirn];
	st_z_rhs[index]     += (st_z_flux[index]     - st_z_flux[indexp1]    ) * dxi[flux_dirn];
	Yet_rhs[index]      += (Yet_flux[index]      - Yet_flux[indexp1]     ) * dxi[flux_dirn];
      }
}
