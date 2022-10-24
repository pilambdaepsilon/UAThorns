#ifndef _RNSID_H_
#define _RNSID_H_ 1

#ifdef __cplusplus
extern "C" 
{
#endif

#ifndef CCTK_REAL
#define CCTK_REAL double
#define CCTK_INT int
#endif
#ifdef RNS_SEQ_COMPILATION
#define cGH int 
#endif
/*-- ifdef cGH -- */
void Hydro_rnsid(const cGH *cctkGH,
             CCTK_REAL *x_grid,
             CCTK_REAL *y_grid,
             CCTK_REAL *z_grid,
             CCTK_REAL eos_k,
             CCTK_REAL eos_ideal_fluid_gamma,
             CCTK_REAL atm_factor,
             CCTK_REAL *Omega_pt,
             CCTK_REAL *R_e_pt,
             CCTK_REAL *r_e_pt,
             CCTK_REAL *mass0_pt,
             CCTK_REAL *gamma_center);
/*-- #endif */

void rnsid(          
          CCTK_REAL rho0_center,
          CCTK_REAL r_ratio,
          CCTK_REAL eos_k,
          CCTK_REAL Gamma_P,
          CCTK_REAL atm_factor,
          CCTK_REAL accuracy,
          char      zero_shift[20],
          char      save_2Dmodel[20],
          char      model2D_file[100],
          char      recover_2Dmodel[20],
          char      rotation_type[20],
          CCTK_REAL A_diff,
	  char eos_type[80],
	  char eos_file[80],
          CCTK_REAL *Omega_pt,
          CCTK_REAL *R_e_pt,
          CCTK_REAL *r_e_pt,
          CCTK_REAL *mass0_pt,
          CCTK_REAL cf,
          CCTK_INT RNS_lmax);

  
  
#ifdef __cplusplus
}
#endif

#endif /* _RNSID_H_ */
