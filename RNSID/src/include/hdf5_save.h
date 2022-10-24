
void hdf5_save_var(int *sdiv, int *mdiv,
		   char fname[],
		   /* ------------------- */
                   char eos_type[],
                   char eos_file[],
                   double *eos_k,
                   double *Gamma_P,
		   /* ------------------- */
                   char   rotation_type[],
                   double *A_diff,
		   /* ------------------- */
                   double *r_ratio,
                   double *rho0_center,	
                   double *r_e,	
		   /* ------------------- */
                   double *s_qp,  /*SDIV+1*/
                   double *mu,	  /*MDIV+1*/
		   /* ------------------- */
                   double **rho, 
                   double **gama, 
                   double **alpha,
                   double **omega,
                   double **energy,
                   double **pressure,
                   double **enthalpy,
                   double **velocity_sq,
                   double *Omega,
                   double *Omega_e,
		   double **Omega_diff);

void hdf5_read_var(int *sdiv, int *mdiv,
		   char fname[],
		   /* ------------------- */
                   char eos_type[],
                   char eos_file[],
                   double *eos_k,
                   double *Gamma_P,
		   /* ------------------- */
                   char   rotation_type[],
                   double *A_diff,
		   /* ------------------- */
                   double *r_ratio,
                   double *rho0_center,
                   double *r_e,	
		   /* ------------------- */
                   double *s_qp,  /*SDIV+1*/
                   double *mu,	  /*MDIV+1*/
		   /* ------------------- */
                   double **rho, 
                   double **gama, 
                   double **alpha,
                   double **omega,
                   double **energy,
                   double **pressure,
                   double **enthalpy,
                   double **velocity_sq,
                   double *Omega,
                   double *Omega_e,
		   double **Omega_diff);

