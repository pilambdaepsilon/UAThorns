#ifndef _RNSID_UTILS_H_
#define _RNSID_UTILS_H_

#define IMAX(a,b) ( a>b ? a : b ) 
#define IMIN(a,b) ( a<b ? a : b ) 
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double **array_allocate(long nrl, long nrh, long ncl, long nch);
void array_free(double **m, long nrl, long nrh, long ncl, long nch);

double ***tensor_allocate(long nrl, 
                   long nrh, 
                   long ncl, 
                   long nch, 
                   long ndl, 
                   long ndh);


void tensor_free(double ***t, 
                   long nrl, 
                   long nrh, 
                   long ncl, 
                   long nch,
               	   long ndl, 
                   long ndh);


double interp_4(double xp[5], 
                double yp[5], 
                int    np ,
                double xb);


void grid_interp(double **old, 
                 double *s_gp, 
                 double *mu, 
                 double r_e, 
                 int nx,
                 int ny,
                 int nz,  
                 double *x_grid,
                 double *y_grid,
                 double *z_grid, 
                 int i,
                 int j,
                 int k, 
                 double *new,
                 int sign);
 

void grid_interp_all( double *s_gp, 
                       double *mu, 
                      double  r_e, 
                      int     nx,
                      int     ny,
                      int     nz, 
                      double *x_grid,
                      double *y_grid,
                      double *z_grid, 
                      int i,
                      int j,
                      int k, 
                      double **nu, 
                      double **B, 
                      double **alpha, 
                      double **omega, 
                      double **nu_dr, 
                      double **B_dr, 
                      double **alpha_dr, 
                      double **omega_dr, 
                      double **nu_dth, 
                      double **B_dth, 
                      double **alpha_dth, 
                      double **omega_dth,
                      double **rho_0, 
                      double **energy, 
                      double **pressure, 
                      double *nu_c, 
                      double *B_c, 
                      double *alpha_c, 
                      double *omega_c, 
                      double *nu_dr_c, 
                      double *B_dr_c, 
                      double *alpha_dr_c, 
                      double *omega_dr_c,                        
                      double *nu_dth_c, 
                      double *B_dth_c, 
                      double *alpha_dth_c, 
                      double *omega_dth_c,
                      double *rho_0_c, 
                      double *energy_c, 
                      double *pressure_c,
                      double *distance_c,
		      double **Omega_diff,
		      double *Omega_diff_c);



void print_arrays_check(
                      int     nx, 
                      int     ny, 
                      double *x_grid, 
                      double *y_grid,
                      int     z_print,
                      double nu_c, 
                      double B_c, 
                      double alpha_c, 
                      double omega_c, 
                      double nu_dr_c, 
                      double B_dr_c, 
                      double alpha_dr_c, 
                      double omega_dr_c,                        
                      double nu_dth_c, 
                      double B_dth_c, 
                      double alpha_dth_c, 
                      double omega_dth_c,
                      double rho_0_c, 
                      double energy_c, 
                      double pressure_c);


     /*
void print_id(
                      int     rdiv, 
                      int     thdiv, 
                      double *rc_grid, 
                      double *thc_grid,
                      double **nu_c, 
                      double **B_c, 
                      double **alpha_c, 
                      double **omega_c, 
                      double **nu_dr, 
                      double **B_dr, 
                      double **alpha_dr, 
                      double **omega_dr,                        
                      double **nu_dth, 
                      double **B_dth, 
                      double **alpha_dth, 
                      double **omega_dth,
                      double **rho_0_c, 
                      double **e_int_c, 
                      double **ut_c, 
                      double **uphi_c);
     */

void transform_units( 
                 char   eos_type[],
                 double  n_P,
                 double  eos_k,
                 double *rho0_center,
                 double *e_center,
                 double *p_center,
                 double *r_e,
                 double **omega,
                 double **energy,
                 double **pressure,
                 double *Mass, 
                 double *Mass_0, 
                 double *T, 
                 double *W, 
                 double *Omega,
                 double *Omega_K,
                 double *R_e,
		 double *Omega_e,
		 double **Omega_diff,
		 double *J);


#endif
