#define _GNU_SOURCE
#include <stdio.h>
#include <unistd.h>
#include "rnsid.h"

#define H5_USE_18_API 
#define H5Acreate_vers 2
#define H5Dcreate_vers 2
#define H5Dopen_vers 2
#include "hdf5.h"
#include "hdf5_save.h"

/* -------
#include <stdio.h>
#include "rnsid.h"
------------ */

#include <string.h> 
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "consts.h"
#include "equil.h"
#include "equil_util.h"
#include "rnsid_util.h"

#ifndef CCTK_REAL
#define CCTK_REAL double
#define CCTK_INT int
#endif

#define AMP " & "

#include "hdf5_save.h"

int main (int argc, char **argv) {
  
  /* ******************************************
     Default Values 
     ********************************************* */
  
  double Omega_pt[1], R_e_pt[1],r_e_pt[1],mass0_pt[1];
  double rho0_center       = 1.28e-3;                /* Central density */
  char   eos_type[80]      = "poly";                 /* {"poly" "tab"} */
  char   eos_file[80]      = "no EOS file specified";
  double eos_k             = 100.0;                  /* poly EOS K */
  double Gamma_P           = 2.0;                    /* poly EOS Gamma */
  char   rotation_type[20] = "diff";                 /* {"uniform" "diff"} */
  double A_diff            = 1.0;                    /* \hat{A} diff rot law */
  double r_ratio           = 1.00;                   /* axes ratio r_p/r_e */
  double accuracy          = 1.e-7;                  /* accuracy goal */
  double atm_factor        = 1.e-10;
  
  char zero_shift[20]      = "no" ;
  char save_2Dmodel[20]    = "no" ;
  char model2D_file_IN[100]  = "INfile.h5";
  char model_file_OUT[100] = "OUTfile.h5";
  char recover_2Dmodel[20] = "no";
  CCTK_REAL cf = 1.0;
  CCTK_INT RNS_lmax = 10;

 /* =============================================== */
  /* Value of the axes ratio used turing iterations  */ 
  /* =============================================== */

  double max_a_ratio   = 1.0;
  double min_a_ratio   = 1.0;
  double delta_a_ratio = 0.1;
 
 
 /* BEGIN RNSID */
  
    char format1[256] ="xx1*| " "%6s"   AMP "%7s"   AMP "%5s"   AMP "%5s"   AMP "%5s"   AMP "%5s"   AMP "%5s"   AMP "%5s"   AMP "%5s"   AMP "%6s"   AMP "%6s\n";
    char format2[256] ="xx2*| " "%6.4f" AMP "%7.5f" AMP "%5.2f" AMP "%5.3f" AMP "%5.3f" AMP "%5.3f" AMP "%5.3f" AMP "%5.3f" AMP "%5.3f" AMP "%6.4f" AMP "%6.4f\n";
    char format4[256] ="gg2*|  %5.3f&%6.4f &%5.1f & %5.2f   &%5.1f&%5.1f&%5.1f&%5.1f& %5.1f & %5.1f & %5.1f & %5.1f & %5.1f & %5.1f & %5.1f \n";
    char format5[256] ="mm1*| " "%6.4f"  " "  "%6.4f"  " "  "%6.4f"  " "  "%7.5f"  " "  "%5.2f"  " "  "%5.3f"  " "  "%5.3f"  " "  "%5.3f"  " "  "%5.3f"  " "  "%5.3f"  " "  "%5.3f"  " "  "%5.3f"  " "  "%5.3f"  " "  "%6.4f"  " "  "%6.4f" "\n" ;

    char format6[256] ="mmA*| %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n";

  /* EQUILIBRIUM VARIABLES */
  
  
  int    n_tab,                        /* Number of points in EOS file */
    a_check=0,                    /* if =200, iteration diverges */ 
    print_dif,                    /* if =1, monitor convergence */
    i;                            /* loop counter */ 
  
  double log_e_tab[MAX_NTAB],     /* energy dens./c^2 in tab. EOS */
    log_p_tab[MAX_NTAB],          /* pressure in tabulated EOS */
    log_h_tab[MAX_NTAB],          /* enthalpy in EOS file */
    log_n0_tab[MAX_NTAB],         /* number density in EOS file */  
    e_center,                     /* central energy density */ 
    p_center,                     /* central pressure */ 
    h_center,                     /* central enthalpy */ 
    e_surface,                    /* surface en. density */ 
    p_surface,                    /* surface pressure */
    enthalpy_min,                 /* minimum enthalpy in EOS */
    s_gp[SDIV+1],                 /* s grid points */
    mu[MDIV+1],                   /* \mu grid points */
    **rho_potential,                          /* potential \rho */ 
    **gama,                         /* potential \gamma */ 
    **omega,                        /* potential \omega */ 
    **alpha,                        /* potential \alpha */ 
    **energy,                       /* energy density \epsilon */
    **pressure,                     /* pressure */ 
    **enthalpy,                     /* enthalpy */
    **velocity_sq,                  /* square of velocity */         
    R_e,                          /* Circumferential radius */
    Mass,                         /* Gravitational mass */
    Mass_0,                       /* Baryon Mass */
    T,                            /* Rotational kinetic energy */
    W,                            /* Gravitational binding energy */
    Omega,                        /* Angular velocity */
    Omega_K,                      /* Ang. vel. of part. in orbit at eq.*/
    r_e,                          /* coord. radius at equator */
    Omega_e,                      /* Ang. vel. at equator, when difrot. */
    **Omega_diff,                   /* Diff. ang. vel. */ 
    J;                            /* Angular momentun */  
  
  /* INITIAL DATA VARIABLES */
  
  
  int m,                             /* counter */
    s,                               /* counter */
    j,                               /* counter */
    k,                               /* counter */
    z_print;                         /* z where check is printed */
  

  double 
   n_P,                           /* polytropic index */
    **nu,                          /* potential nu */
    **B,                           /* potential B */  
    **rho_0,                       /* rest mass density */
    **nu_dr,                       /* r-der. in s-coord. of nu */
    **B_dr,                        /* r-der. in s-coord. of B */
    **alpha_dr,                    /* r-der. in s-coord. of alpha */
    **omega_dr,                    /* r-der. in s-coord. of omega */
    **nu_dth,                      /* theta-der. in mu-coord. of nu */
    **B_dth,                       /* theta-der. in mu-coord. of B */
    **alpha_dth,                   /* theta-der. in mu-coord. of alpha */
    **omega_dth,                   /* theta-der. in mu-coord. of omega */
    x_i,                           /* x at i */
    y_j,                           /* y at j */
    z_k,                           /* z at k */       
    nu_ijk,                        /* nu at ijk point */
    exp_nu_ijk,                    /* exp(nu) at ijk point */
    B_ijk,                         /* B at ijk point */
    omega_ijk,                     /* omega at ijk point */
    alpha_ijk,                     /* alpha at ijk point */
    exp_alpha_ijk,                 /* exp(alpha) at ijk point */
    rho_0_ijk,                     /* rho_0 at ijk point */
    energy_ijk,                    /* energy at ijk point */
    pressure_ijk,                  /* pressure at ijk point */
    nu_dx,                         /* derivative of nu w.r.t. x */
    nu_dy,                         /* derivative of nu w.r.t. y */
    B_dx,                          /* derivative of B w.r.t. x */
    B_dy,                          /* derivative of B w.r.t. y */
    omega_dx,                      /* derivative of omega w.r.t. x */
    omega_dy,                      /* derivative of omega w.r.t. y */
    omega_dz,                      /* derivative of omega w.r.t. z */
    alpha_dx,                      /* derivative of alpha w.r.t. x */
    alpha_dy,                      /* derivative of alpha w.r.t. y */
    r_ijk,                         /* r at ijk point */
    r_bar_ijk,                     /* sqrt(x^2+y^2) at ijk point */
    dr_dx,                         /* dr/dx */
    dr_dy,                         /* dr/dy */
    dr_dz,                         /* dr/dz */
    dtheta_dx,                     /* dtheta/dx */
    dtheta_dy,                     /* dtheta/dy */
    dtheta_dz,                     /* dtheta/dz */
    nu_dr_ijk,                     /* dnu/dr at ijk */
    B_dr_ijk,                      /* dB/dr at ijk */
    alpha_dr_ijk,                  /* dalpha/dr at ijk */
    omega_dr_ijk,                  /* domega/dr at ijk */
    nu_dtheta_ijk,                 /* dnu/dtheta at ijk */
    B_dtheta_ijk,                  /* dB/dtheta at ijk */
    alpha_dtheta_ijk,              /* dalpha/dtheta at ijk */
    omega_dtheta_ijk,              /* domega/dtheta at ijk */
    gamma_ijk,                     /* gamma = det(3g) */
    W_ijk,                         /* Lorentz factor */     
    h_ijk,                         /* h = 1 + eps + P/rho */
    distance_ijk,                  /* Signed distance to surface */
    e_atm,                         /* energy density of atmosphere */
    p_atm,                         /* pressure of atmosphere */
    rho_0_atm,                     /* rest mass density of atmosphere */
    dens_atm,                      /* D of atmosphere */
    tau_atm,                       /* tau of atmosphere */
    temp_a,                        /* temporary variables */
    temp_o,
    temp_g,
    temp_r,
    temp_e,
    temp_p,
    temp_h,
    temp_v,
    Omega_ijk;

  int n_nearest;
  double n0;
  int ii,jj,kk,idx;
  /* =============================================== */
  /* Intermediate quantities used for output purpose */ 
  /* =============================================== */
  double Pa,Pe,td,beta,M_over_R;
  double invPa,invPe,invtd;
  /// Total number of grid points
  int np ;
  // -----------------------------
  //  GRID size  (default values)
  // -----------------------------
  int     Px,Py,Pz,Nx,Ny,Nz;
  double  Lx,Ly,Lz,dx;
  Px = Py = Pz = 25; 
  Nx = Ny = Nz = 49;
  Lx = Ly = Lz = 12;
  dx = 0.5;
  np = Nx*Ny*Nz;


  /* ******************************************
     .... a la charte
     ********************************************* */
  int opterr,c;
  opterr = 0;
  for(c = getopt(argc,argv,"+ht:f:k:g:i:o:x:y:z:d:");
      c != EOF;
      c = getopt(argc,argv,"+ht:f:k:g:i:o:x:y:z:d:")) {
    switch(c) {
    case 'x' : /* Number of grid point x-direction from 0 (included) to BORDER */
      sscanf(optarg,"%d",&Px);
      break;
    case 'y' : /* Number of grid point x-direction from 0 (included) to BORDER */
      sscanf(optarg,"%d",&Py);
      break;
    case 'z' : /* Number of grid point x-direction from 0 (included) to BORDER */
      sscanf(optarg,"%d",&Pz);
      break;
    case 'd' : /* grid spacing in units G=c=M_sun=1 */
      sscanf(optarg,"%lf",&dx);
      break;
    case 'i' : /* initial_model_data.h5  */
      sscanf(optarg,"%s",model2D_file_IN);
      break;
    case 'o' : /* ginal_model_data.h5 */
      sscanf(optarg,"%s",model_file_OUT);
      break;
    case 't' : /* eos type */
      sscanf(optarg,"%s",eos_type);
      break;
    case 'f' : /* eos TAB file */
      sscanf(optarg,"%s",eos_file);
      break;
    case 'k' : /* poly EOS K */
      sscanf(optarg,"%lf",&eos_k);
      break;
    case 'g': /* poly EOS Gamma */
      sscanf(optarg,"%lf",&Gamma_P);
      break;
    case 'h':
      printf("\n");
      printf("Help: ./RNS_readID -options<meaning> {default value} \n");
      printf("\n");
      printf("      -i <initial_model_data.hdf> \n");
      printf("      -o <initial_model_grid_data.hdf> \n");
      printf("\n");
      printf("      -d <step-in-coordinates> \n");
      printf("      -x <number of points in x-dir>   \n");
      printf("      -y <number of points in x-dir>   \n");
      printf("      -z <number of points in x-dir>   \n");
      printf("      [GRID size (d [-Nz,Nz]) x (d [-Ny,Ny])x (d [-Nx,Nx])] \n"); 
      printf("\n");
      printf("      -t <EOS type: poly/tab>  {poly} \n");
      printf("      -f <tab EOS file>        {no file} \n");
      printf("      -k <poly EOS K>          {100.0}\n");
      printf("      -g <poly EOS Gamma>      {2.0}\n");
      printf("\n");
      return 0;
    default:
      break;
    }
  }


  Nx = 2*Px -1; 
  Ny = 2*Py -1; 
  Nz = 2*Pz -1; 
  Lx = (Px-1) * dx;
  Ly = (Py-1) * dx;
  Lz = (Pz-1) * dx;
  np =  Nx*Ny*Nz;
  printf("#====================================================\n");
  printf("# Grid setting x range is: [ %e , %e ]\n",-Lx,Lx);
  printf("#              y range is: [ %e , %e ]\n",-Ly,Ly);
  printf("#              z range is: [ %e , %e ]\n",-Lz,Lz);
  printf("# Number of points is (%d) %d,%d,%d \n",np,Nz,Ny,Nx) ;
  printf("#====================================================\n");

  /* ******************************************
     call rnsid
     ********************************************* */

  /* Read valus from tabulated EOS and convert to C.U. */  

  if(strcmp(eos_type,"tab")==0) {
    printf(" TAB eos from file: %s\n",eos_file);    
    load_eos( eos_file, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, &n_tab );
  }  
 
  /* COMPUTE POLYTROPIC INDEX AND CENTRAL ENERGY DENSITY */
  
  n_P=1.0/(Gamma_P-1.0);            

  if(strcmp(eos_type,"tab")==0) {
      /*this should work now*/
      n_nearest = 50;
      n0 = rho0_center/(MB*cactusM);
      e_center = pow(10.0,interp(log_n0_tab, log_e_tab, n_tab,log10(n0), &n_nearest));
  }
  else {
      e_center = (eos_k*pow(rho0_center,Gamma_P)/(Gamma_P-1.0)+rho0_center);
  }
 

  /* SET UP GRID */
  make_grid(s_gp, mu);
  

  /* ALLLOCATE MEMORY */
  rho_potential = array_allocate(1,SDIV,1,MDIV);
  gama = array_allocate(1,SDIV,1,MDIV);
  alpha = array_allocate(1,SDIV,1,MDIV);
  omega = array_allocate(1,SDIV,1,MDIV);
  energy = array_allocate(1,SDIV,1,MDIV);
  pressure = array_allocate(1,SDIV,1,MDIV);
  enthalpy = array_allocate(1,SDIV,1,MDIV);
  velocity_sq = array_allocate(1,SDIV,1,MDIV);
  Omega_diff = array_allocate(1,SDIV,1,MDIV);
  

  /* SET DEFAULT EQUILIBRIUM PARAMETERS */
  if(strcmp(eos_type,"tab")==0) {
      /*to use C.U. here different values have to be set*/
    e_surface=7.8e-15;
    p_surface=1.12379e-28;
    enthalpy_min=1.0/(C*C);
    printf("e_surface = %g\n",e_surface);
    printf("p_surface = %g\n",p_surface);
    
  }
  else {
    e_surface=0.0;
    p_surface=0.0;
    enthalpy_min=0.0;
  }
  
  Omega_e=0.0; /* initialize ang. vel. at equator for diff. rot. */
  
  print_dif=1;
  
  int sdiv=SDIV;
  int mdiv=MDIV;
  printf("Reading data from file: %s\n",model2D_file_IN);
  hdf5_read_var(&sdiv,&mdiv,model2D_file_IN,eos_type,eos_file,&eos_k,&Gamma_P,rotation_type,&A_diff,&r_ratio,&rho0_center,&r_e,
                s_gp, mu,rho_potential, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq, &Omega,&Omega_e,Omega_diff);
  comp_values( s_gp, mu, r_ratio, e_surface, r_e, eos_type, log_e_tab,
  		   log_n0_tab, n_tab, Omega, rho_potential, gama, alpha, omega, 
  	       energy, pressure, enthalpy, velocity_sq, &Mass, 
  	       &Mass_0, &T, &W, &Omega_K, &R_e, rotation_type, Omega_diff,
  	       &J);

  invPa = Omega_diff[1][1]         * (cactusT/(1000*2*PI));
  invPe = Omega_e                  * (cactusT/(1000*2*PI));
  invtd = sqrt(R_e*R_e*R_e/Mass_0) * (cactusT/(1000*2*PI));
  beta     = T/W;
  M_over_R = Mass/R_e;
  printf("xx0*|============================================================================================\n");	
  printf(format1,"rho_c","r_p/r_e","R_e","M_0","M","J","P_a" ,"P_e" ,"t_D"  ,"M/R_e","beta");
  printf(format1,"(e-4)",""       ,""   ,""   ,"" ,"" ,"msec","msec","msec" ,""     ,"");
  printf("xx0*|--------------------------------------------------------------------------------------------\n");	
  printf(format2,10000.0*rho0_center,r_ratio,R_e,Mass_0,Mass,J,1.0/invPa,1.0/invPe,1.0/invtd,M_over_R,beta);
  printf("xx3*|============================================================================================\n");
  


  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************
  /* CONSTRUCT ARRAYS WITH NEEDED POLAR QUANTITIES */
  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************

  /// 1-D array storing the values of coordinate x of the {\tt np} grid points [unit: km]
  double* x_grid ;
  /// 1-D array storing the values of coordinate y of the {\tt np} grid points [unit: km]
  double* y_grid ;
  /// 1-D array storing the values of coordinate z of the {\tt np} grid points [unit: km]
  double* z_grid ;
  /// Lapse function $N$ at the {\tt np} grid points (1-D array)
  double* alp   ;   /// Component $\beta^x$ of the shift vector of non rotating coordinates [unit: $c$]
  double* betax ; /// Component $\beta^y$ of the shift vector of non rotating coordinates [unit: $c$]
  double* betay ; /// Component $\beta^z$ of the shift vector of non rotating coordinates [unit: $c$]
  double* betaz ; /// Metric coefficient $\gamma_{xx}$ at the grid points (1-D array)
  double* gxx ; 
  /// Metric coefficient $\gamma_{xy}$ at the grid points (1-D array)
  double* gxy ;
  /// Metric coefficient $\gamma_{xz}$ at the grid points (1-D array)
  double* gxz ; 
  /// Metric coefficient $\gamma_{yy}$ at the grid points (1-D array)
  double* gyy ; 
  /// Metric coefficient $\gamma_{yz}$ at the grid points (1-D array)
  double* gyz ; 
  /// Metric coefficient $\gamma_{zz}$ at the grid points (1-D array)
  double* gzz ; 
  /// Component $K_{xx}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
  double* kxx ;
  /// Component $K_{xy}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
  double* kxy ;
  /// Component $K_{xz}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
  double* kxz ;
  /// Component $K_{yy}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
  double* kyy ;
  /// Component $K_{yz}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
  double* kyz ;
  /// Component $K_{zz}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
  double* kzz ;
  // Hydro components
  //------------------
  /** Baryon density in the fluid frame at the {\tt np} grid points (1-D array)
   * [unit: ${\rm kg \, m}^{-3}$]
   */
  double* rho ;
  double* press ;
  /// Specific internal energy at the  {\tt np} grid points (1-D array) [unit: $c^2$]
  double* eps ;
  /** Component $U^x$ of the fluid 3-velocity with respect to the Eulerian
   * observer, at the {\tt np} grid points (1-D array) [unit: $c$]
   */
  double* velx ;
  /** Component $U^y$ of the fluid 3-velocity with respect to the Eulerian
   * observer, at the {\tt np} grid points (1-D array) [unit: $c$]
   */
  double* vely ;
  /** Component $U^z$ of the fluid 3-velocity with respect to the Eulerian
   * observer, at the {\tt np} grid points (1-D array) [unit: $c$]
   */
  double* velz ;


  x_grid = malloc(np*sizeof(double)); // new double[np] ;
  y_grid = malloc(np*sizeof(double)); // new double[np] ;
  z_grid = malloc(np*sizeof(double)); // new double[np] ;

  alp = malloc(np*sizeof(double)); // new double[np] ;

  betax = malloc(np*sizeof(double)); // new double[np] ;
  betay = malloc(np*sizeof(double)); // new double[np] ;
  betaz = malloc(np*sizeof(double)); // new double[np] ;

  gxx = malloc(np*sizeof(double)); // new double[np] ;
  gxy = malloc(np*sizeof(double)); // new double[np] ;
  gxz = malloc(np*sizeof(double)); // new double[np] ;
  gyy = malloc(np*sizeof(double)); // new double[np] ;
  gyz = malloc(np*sizeof(double)); // new double[np] ;
  gzz = malloc(np*sizeof(double)); // new double[np] ;

  kxx = malloc(np*sizeof(double)); // new double[np] ;
  kxy = malloc(np*sizeof(double)); // new double[np] ;
  kxz = malloc(np*sizeof(double)); // new double[np] ;
  kyy = malloc(np*sizeof(double)); // new double[np] ;
  kyz = malloc(np*sizeof(double)); // new double[np] ;
  kzz = malloc(np*sizeof(double)); // new double[np] ;

  rho   = malloc(np*sizeof(double)); // new double[np] ;
  press = malloc(np*sizeof(double)); // new double[np] ;
  eps   = malloc(np*sizeof(double)); // new double[np] ;

  velx = malloc(np*sizeof(double)); // new double[np] ;
  vely = malloc(np*sizeof(double)); // new double[np] ;
  velz = malloc(np*sizeof(double)); // new double[np] ;


  for (ii=0; ii < Nx ; ii++) 
    for (jj=0; jj < Ny ; jj++) 
      for (kk=0; kk < Nz ; kk++) {
	idx = ii + jj * Nx + kk * Nx*Ny;
        x_grid[idx] = -Lx + dx * ii;
        y_grid[idx] = -Ly + dx * jj;
        z_grid[idx] = -Lz + dx * kk;
        // x_grid[idx] = -Lx + (2*Lx/(Nx-1)) * i;
        // y_grid[idx] = -Ly + (2*Ly/(Ny-1)) * j;
        // z_grid[idx] = -Lz + (2*Lz/(Nz-1)) * k;
      } 

  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************
  //**************************************************************************************************

  nu    = array_allocate(1,SDIV,1,MDIV);
  B     = array_allocate(1,SDIV,1,MDIV);
  rho_0 = array_allocate(1,SDIV,1,MDIV);      

  for(m=1;m<=MDIV;m++) 
     for(s=1;s<=SDIV;s++) {
  	      nu[s][m] = (gama[s][m]+rho_potential[s][m])/2.0; 
  	      B[s][m] = exp(gama[s][m]);	       
          rho_0[s][m] = (energy[s][m]+pressure[s][m])
                         *exp(-enthalpy[s][m]);
     } 
  
  array_free(rho_potential,1,SDIV,1,MDIV);
  array_free(gama,1,SDIV,1,MDIV);
  array_free(enthalpy,1,SDIV,1,MDIV);
  array_free(velocity_sq,1,SDIV,1,MDIV);
  
  nu_dr = array_allocate(1,SDIV,1,MDIV);
  B_dr = array_allocate(1,SDIV,1,MDIV);
  alpha_dr = array_allocate(1,SDIV,1,MDIV);
  omega_dr = array_allocate(1,SDIV,1,MDIV);
  nu_dth = array_allocate(1,SDIV,1,MDIV);
  B_dth = array_allocate(1,SDIV,1,MDIV);
  alpha_dth = array_allocate(1,SDIV,1,MDIV);
  omega_dth = array_allocate(1,SDIV,1,MDIV);

  for(m=1;m<=MDIV;m++) 
     for(s=1;s<=SDIV;s++) {
        nu_dr[s][m] = deriv_s(nu,s,m)*SQ(1.0-s_gp[s])/r_e;
        B_dr[s][m] = deriv_s(B,s,m)*SQ(1.0-s_gp[s])/r_e;
        alpha_dr[s][m] = deriv_s(alpha,s,m)*SQ(1.0-s_gp[s])/r_e;
        omega_dr[s][m] = deriv_s(omega,s,m)*SQ(1.0-s_gp[s])/r_e;
        nu_dth[s][m] = deriv_m(nu,s,m)*(-sqrt(1.0-SQ(mu[m])));
        B_dth[s][m] = deriv_m(B,s,m)*(-sqrt(1.0-SQ(mu[m])));
        alpha_dth[s][m] = deriv_m(alpha,s,m)*(-sqrt(1.0-SQ(mu[m])));
        omega_dth[s][m] = deriv_m(omega,s,m)*(-sqrt(1.0-SQ(mu[m])));
     }  

  /* COMPUTE INITIAL DATA */

  // rho_0_atm = rnsid_rho_min; /* rename the constant for historical reasons */
  rho_0_atm = 1e-10; /* rename the constant for historical reasons */
  e_atm = rho_0_atm;
  p_atm = eos_k*pow(rho_0_atm,Gamma_P);

  
  for(i=1;i<=Nx;i++)
     for(j=1;j<=Ny;j++) 
       for(k=1;k<=Nz;k++) {

           x_i = x_grid[i-1+Nx*(j-1+Ny*(k-1))];
           y_j = y_grid[i-1+Nx*(j-1+Ny*(k-1))];
           z_k = z_grid[i-1+Nx*(j-1+Ny*(k-1))];

           grid_interp_all( s_gp, mu, r_e, Nx, Ny, Nz, 
                            x_grid, y_grid, z_grid,
                            i, j, k,
                            nu, B, alpha, omega, 
                            nu_dr, B_dr, alpha_dr, omega_dr, 
                            nu_dth, B_dth, alpha_dth, omega_dth,
                            rho_0, energy, pressure,
                            &nu_ijk, &B_ijk, &alpha_ijk, &omega_ijk, 
                            &nu_dr_ijk, &B_dr_ijk, &alpha_dr_ijk, 
                            &omega_dr_ijk,             
                            &nu_dtheta_ijk, &B_dtheta_ijk, 
                            &alpha_dtheta_ijk, &omega_dtheta_ijk,
                            &rho_0_ijk, &energy_ijk, &pressure_ijk,
                            &distance_ijk,
  			    Omega_diff, &Omega_ijk);

          /* *************************************** */
  	      /* DETECT if it is in the ATMOSPHERE       */
          /* *************************************** */
          if( (rho_0_ijk<=0.0) || (energy_ijk<=0.0) || 
              (pressure_ijk<=0.0) ) {
 
            rho_0_ijk   = rho_0_atm;
            energy_ijk  = e_atm + 1.e-20;
            pressure_ijk= p_atm;
            
          }
          /* *************************************** */
  	      /* END ATMOSPERE SETTINGs                      */
          /* *************************************** */
  	      
          exp_nu_ijk = exp(nu_ijk);
          exp_alpha_ijk = exp(alpha_ijk);
 
          r_ijk = sqrt(SQ(x_i)+SQ(y_j)+SQ(z_k));
          r_bar_ijk = sqrt(SQ(x_i)+SQ(y_j));


          alp[i-1+Nx*(j-1+Ny*(k-1))] = exp_nu_ijk;


          if(x_i==0.0 && y_j==0.0) {
  
            gxx[i-1+Nx*(j-1+Ny*(k-1))] = SQ(exp_alpha_ijk);
            gyy[i-1+Nx*(j-1+Ny*(k-1))] = SQ(exp_alpha_ijk);
            gzz[i-1+Nx*(j-1+Ny*(k-1))] = SQ(exp_alpha_ijk);

            gxy[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
            gxz[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
            gyz[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
  
            kxx[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
            kyy[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
            kzz[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
            kxy[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
            kxz[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
            kyz[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;

          } else {
 
            dr_dx = x_i / r_ijk;
            dr_dy = y_j / r_ijk;
            dr_dz = z_k / r_ijk;

            dtheta_dx = x_i*z_k/(SQ(r_ijk)*r_bar_ijk);
            dtheta_dy = y_j*z_k/(SQ(r_ijk)*r_bar_ijk);
            dtheta_dz = - r_bar_ijk/SQ(r_ijk);
 
            nu_dx = dr_dx*nu_dr_ijk + dtheta_dx*nu_dtheta_ijk;
            nu_dy = dr_dy*nu_dr_ijk + dtheta_dy*nu_dtheta_ijk;

            B_dx = dr_dx*B_dr_ijk + dtheta_dx*B_dtheta_ijk;
            B_dy = dr_dy*B_dr_ijk + dtheta_dy*B_dtheta_ijk;
 
            alpha_dx = dr_dx*alpha_dr_ijk + dtheta_dx*alpha_dtheta_ijk;
            alpha_dy = dr_dy*alpha_dr_ijk + dtheta_dy*alpha_dtheta_ijk;
  
            omega_dx = dr_dx*omega_dr_ijk + dtheta_dx*omega_dtheta_ijk;
            omega_dy = dr_dy*omega_dr_ijk + dtheta_dy*omega_dtheta_ijk;

  		/* enforce omega_dz=0 at z=0 (it is slightly noNzero due
  		   to O(h) forwards formula in computing derivative) */

            if(z_k==0.0) 
              omega_dz = 0.0;
            else
              omega_dz = dr_dz*omega_dr_ijk + dtheta_dz*omega_dtheta_ijk;

 
            gxx[i-1+Nx*(j-1+Ny*(k-1))] = ( SQ(B_ijk*y_j/exp_nu_ijk) 
                                         +SQ(exp_alpha_ijk*x_i) ) /
                                         (SQ(x_i)+SQ(y_j));

            gxy[i-1+Nx*(j-1+Ny*(k-1))] = ( SQ(exp_alpha_ijk) 
                                        -SQ(B_ijk/exp_nu_ijk) ) *
                                        x_i*y_j/(SQ(x_i)+SQ(y_j));

            gxz[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;

            gyy[i-1+Nx*(j-1+Ny*(k-1))] = ( SQ(B_ijk*x_i/exp_nu_ijk) 
                                         +SQ(exp_alpha_ijk*y_j) ) /
                                         (SQ(x_i)+SQ(y_j));
  
            gyz[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
  
            gzz[i-1+Nx*(j-1+Ny*(k-1))] = SQ(exp_alpha_ijk);
 

            kxx[i-1+Nx*(j-1+Ny*(k-1))] = (  ( SQ(r_bar_ijk)*y_j*omega_dx+
                                         (x_i*nu_dy-y_j*nu_dx)*SQ(y_j)
                                         *omega_ijk)*SQ(B_ijk)         
                                         +(y_j*B_dx-x_i*B_dy)*omega_ijk
                                         *SQ(y_j)*B_ijk
                                         +(y_j*alpha_dx-x_i*alpha_dy)
                                         *omega_ijk*SQ(x_i*exp_alpha_ijk
                                         *exp_nu_ijk))/(SQ(r_bar_ijk
                                         *exp_nu_ijk)*exp_nu_ijk);
  		
            kxy[i-1+Nx*(j-1+Ny*(k-1))] = ( ( 0.5*SQ(r_bar_ijk)*
                                         (y_j*omega_dy - x_i*omega_dx) + 
                                         (y_j*nu_dx-x_i*nu_dy)*x_i*y_j*
                                         omega_ijk )*SQ(B_ijk) 
                                         +(-y_j*B_dx+x_i*B_dy)*omega_ijk
                                         *x_i*y_j*B_ijk
                                         +(y_j*alpha_dx-x_i*alpha_dy)
                                         *omega_ijk*x_i*y_j
                                         *SQ(exp_alpha_ijk*exp_nu_ijk))/
                                         (SQ(r_bar_ijk*exp_nu_ijk)
                                         *exp_nu_ijk);
  	        
            kxz[i-1+Nx*(j-1+Ny*(k-1))] = 0.5*SQ(B_ijk)*y_j*omega_dz/
                                         ( SQ(exp_nu_ijk)*exp_nu_ijk );
  
            kyy[i-1+Nx*(j-1+Ny*(k-1))] = ( ( -SQ(r_bar_ijk)*x_i*omega_dy+
                                         (x_i*nu_dy-y_j*nu_dx)*SQ(x_i)* 
                                         omega_ijk )*SQ(B_ijk) 
                                         +(y_j*B_dx-x_i*B_dy)*omega_ijk
                                         *SQ(x_i)*B_ijk
                                         +(y_j*alpha_dx-x_i*alpha_dy)
                                         *omega_ijk*SQ(y_j*exp_alpha_ijk
                                         *exp_nu_ijk))/(SQ(r_bar_ijk
                                         *exp_nu_ijk)*exp_nu_ijk);
  	       
            kyz[i-1+Nx*(j-1+Ny*(k-1))] = -0.5*SQ(B_ijk)*x_i*omega_dz/
                                         ( SQ(exp_nu_ijk)*exp_nu_ijk );
 
            kzz[i-1+Nx*(j-1+Ny*(k-1))] = (y_j*alpha_dx-x_i*alpha_dy)*
                                         omega_ijk*SQ(exp_alpha_ijk)/
                                         exp_nu_ijk;
  	      }
 
  	       
          betax[i-1+Nx*(j-1+Ny*(k-1))] = omega_ijk*y_j;
          betay[i-1+Nx*(j-1+Ny*(k-1))] = -omega_ijk*x_i;
  	       
          betaz[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
  	       
          rho[i-1+Nx*(j-1+Ny*(k-1))] = rho_0_ijk;
  	       
          eps[i-1+Nx*(j-1+Ny*(k-1))] = energy_ijk/rho_0_ijk-1.0;
          h_ijk = (energy_ijk+pressure_ijk)/rho_0_ijk;
         
          gamma_ijk = SQ(exp_alpha_ijk)*SQ(B_ijk*exp_alpha_ijk/
                      exp_nu_ijk); 

          W_ijk = 1.0/sqrt(1.0-SQ((omega_ijk-Omega_ijk)*B_ijk*
                  r_bar_ijk/SQ(exp_nu_ijk)));

//          w_lorentz[i-1+Nx*(j-1+Ny*(k-1))] = W_ijk;
//  
//          dens[i-1+Nx*(j-1+Ny*(k-1))] = sqrt(gamma_ijk)*W_ijk*rho_0_ijk;
//
//          tau[i-1+Nx*(j-1+Ny*(k-1))] = sqrt(gamma_ijk)*( rho_0_ijk*h_ijk
//                                       *SQ(W_ijk)-pressure_ijk - 
//                                       W_ijk*rho_0_ijk );
//
//  	      sx[i-1+Nx*(j-1+Ny*(k-1))] = sqrt(gamma_ijk)*rho_0_ijk*h_ijk
//                                      *SQ(W_ijk*B_ijk/exp_nu_ijk)*
//                                      (omega_ijk-Omega_ijk)*y_j/exp_nu_ijk;
//
//          sy[i-1+Nx*(j-1+Ny*(k-1))] = - sqrt(gamma_ijk)*rho_0_ijk*h_ijk
//                                      *SQ(W_ijk*B_ijk/exp_nu_ijk)*
//                                      (omega_ijk-Omega_ijk)*x_i/exp_nu_ijk;
//
//          sz[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
//
  	  velx[i-1+Nx*(j-1+Ny*(k-1))] = (omega_ijk-Omega_ijk)
  	                                     *y_j/exp_nu_ijk;

          vely[i-1+Nx*(j-1+Ny*(k-1))] = -(omega_ijk-Omega_ijk)
  		                              *x_i/exp_nu_ijk;

          velz[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;

          press[i-1+Nx*(j-1+Ny*(k-1))] =pressure_ijk;
  	      
  	  dens_atm = sqrt(gamma_ijk)*rho_0_atm;
          tau_atm = sqrt(gamma_ijk)*eos_k*pow(rho_0_atm,Gamma_P) / 
            (Gamma_P - 1.0);

          /* *************************************** */
  	      /* ATMOSPERE SETTINGs                      */
          /* *************************************** */
	  //if ( (rho[i-1+Nx*(j-1+Ny*(k-1))] < (1.0 + RNS_atmo_tolerance) * rho_0_atm) ) {
  	      if ( (rho[i-1+Nx*(j-1+Ny*(k-1))] < (1.0 + 1e-3) * rho_0_atm) ) {
  		velx[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
  		vely[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
  		velz[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
  		// dens[i-1+Nx*(j-1+Ny*(k-1))] = dens_atm;
  		// tau[i-1+Nx*(j-1+Ny*(k-1))] = tau_atm;
  		// sx[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
  		// sy[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
  		// sz[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
  	      } 
          /* *************************************** */
  	      /* ATMOSPERE SETTINGs                      */
          /* *************************************** */


  	   } /* END FOR LOOP OF LINE 525 */




  //if( strcmp(zero_shift, "yes")==0) {
  //
  //  /* SET SHIFT TO ZERO */
  //
  //  for(i=1;i<=Nx;i++)
  //     for(j=1;j<=Ny;j++)
  //        for(k=1;k<=Nz;k++) {
  //           betax[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
  //           betay[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;     
  //           betaz[i-1+Nx*(j-1+Ny*(k-1))] = 0.0;
  //        }
  //}

  // ------------------------------------------
  // HDF5 save
  // ------------------------------------------

  hid_t       file_id;   /* file identifier */
  herr_t      status;

  hid_t       attribute_id; /* identifiers for attributes*/
  hid_t       attributeH5type,attributeH5c128;

  hid_t       dataset_id, dataspace_id;  /* identifiers for dsets*/
  hsize_t     dims[3];

  //int         ii,jj;
  double      *dset_data;
  double      **var;
  int         varIndex;

  //cout << "I am here " << endl;

  /* =========================================== */
  /* Create a new file using default properties. */
  /* =========================================== */
  file_id = H5Fcreate(model_file_OUT, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  dims[2] = Nx;
  dims[1] = Ny;
  dims[0] = Nz;

  double* varlist[]  = {x_grid,y_grid,z_grid,
                        alp,betax,betay,betaz,
                        gxx,gxy,gxz,gyy,gyz,gzz,
                        kxx,kxy,kxz,kyy,kyz,kzz,
                        rho,press,eps,velx,vely,velz
                       }; 
  char  * varnames[] = {"/X[0]","/X[1]","/X[2]",
                        "/alp","/betax","/betay","/betaz",
			"/gxx","/gxy","/gxz","/gyy","/gyz","/gzz",
			"/Kxx","/Kxy","/Kxz","/Kyy","/Kyz","/Kzz",
                        "/rho","/press","/eps","/vel[0]","/vel[1]","/vel[2]"
                       }; 

  for(idx=0;idx< 25; idx++) {
    dataspace_id = H5Screate_simple(3, dims, NULL);
    dataset_id = H5Dcreate(file_id, varnames[idx], 
    			 H5T_NATIVE_DOUBLE, dataspace_id, 
    			 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
    		    varlist[idx]);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
  }

  status = H5Fclose(file_id); 

  // ------------------------------------------
  // Deallocate variables
  // ------------------------------------------
  /* compute central value of 3-determinant of metric */
  // *gamma_center = SQ(B[1][1])*exp( 4.0*alpha[1][1] - 2.0*nu[1][1]);
  /* FREE MEMORY */
  	      
  array_free(alpha,1,SDIV,1,MDIV);
  array_free(omega,1,SDIV,1,MDIV);      
  array_free(rho_0,1,SDIV,1,MDIV);      
  array_free(energy,1,SDIV,1,MDIV);
  array_free(pressure,1,SDIV,1,MDIV);

  array_free(nu,1,SDIV,1,MDIV);
  array_free(B,1,SDIV,1,MDIV);

  array_free(nu_dr,1,SDIV,1,MDIV);
  array_free(B_dr,1,SDIV,1,MDIV);
  array_free(alpha_dr,1,SDIV,1,MDIV);
  array_free(omega_dr,1,SDIV,1,MDIV);
  array_free(nu_dth,1,SDIV,1,MDIV);
  array_free(B_dth,1,SDIV,1,MDIV);
  array_free(alpha_dth,1,SDIV,1,MDIV);
  array_free(omega_dth,1,SDIV,1,MDIV);

  //free(s_gp);
  //free(mu);
    

  return 0;
}
