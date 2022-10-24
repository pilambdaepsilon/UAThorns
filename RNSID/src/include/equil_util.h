
void hunt(double xx[], int n, double x, int *jlo);

double interp(double xp[], 
              double yp[], 
              int    np ,
              double xb, 
              int    *n_nearest_pt);

double deriv_s(double **f,int s, int m);

double deriv_ss(double **f,int s, int m);

double deriv_m(double **f,int s, int m);

double deriv_mm(double **f,int s, int m);

double deriv_sm(double **f,int s, int m);

double legendre( int n, double x );

double plgndr(int l, int m, double x);
 
double rtsec_G( double (*func)(double, double, double), 
                double Gamma_P,
                double x1, 
                double x2, 
                double xacc,
                double ee,
		double eos_k);

double zbrent_diff(double (*func)(double, double, double, double, double, 
                                  double, double, double), 
                   double r_e,
                   double rho_equator_h,
                   double gama_equator_h,
                   double omega_equator_h,
                   double rho_pole_h,
                   double gama_pole_h,  
                   double x1, 
                   double x2, 
                   double tol,
		   double A_diff);

double zbrent_rot(double (*func)(double, double, double, double, double, 
                                  double, double, double), 
                   double r_e,
                   double rhogp, 
                   double omegagp, 
                   double sgp, 
                   double mugp, 
                   double Omega_c,
                   double x1, 
                   double x2, 
    		   double tol,
		   double A_diff);

double diff_rotation(double x, 
                     double r_e, 
                     double rho_equator_h, 
                     double gama_equator_h, 
                     double omega_equator_h, 
                     double rho_pole_h, 
                     double gama_pole_h,
		     double A_diff);

double rotation_law( double x, 
                     double r_e, 
                     double rhogp, 
                     double omegagp, 
                     double sgp, 
                     double mugp, 
                     double Omega_c,
		     double A_diff);
