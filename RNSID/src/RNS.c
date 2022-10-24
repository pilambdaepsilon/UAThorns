#define _GNU_SOURCE
#include <stdio.h>
#include <unistd.h>
#include "rnsid.h"

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
  char model2D_file_IN[100]  = ""; //"INfile.h5";
  char model2D_file_OUT[100] = ""; //"OUTfile.h5";
  char recover_2Dmodel[20] = "no";
  char export_data_file_name[100] = "";
  FILE * export_data_file;  
  int sdiv=SDIV;
  int mdiv=MDIV;

  CCTK_REAL cf = 1.0;
  CCTK_INT RNS_lmax = 10;


  /* =============================================== */
  /* Value of the axes ratio used during iterations  */ 
  /* =============================================== */
  double max_a_ratio   = 1.0;
  double min_a_ratio   = 1.0;
  double delta_a_ratio = 0.1;
  /* =============================================== */
  /* Value of the rho0 used during iterations        */
  /* =============================================== */
  double max_rho0_center   = 1.28e-3;
  double min_rho0_center   = 1.28e-3;
  double delta_rho0_center = 0.1e-3;
  double rho0_centerITR;
  int    idx = 0;
  int    ii; 
 
  /* =============================================== */
  /* FORMAT STRING USED FOR OUTPUT                   */
  /* =============================================== */

  /* FORMAT FOR A NICE TABLE TO BE IMPORTED IN LATEX (11 itemes)*/ 
  char format1[256] ="xx1*| " "%6s"   AMP "%7s"   AMP "%5s"   AMP "%5s"   AMP "%5s"   AMP "%5s"   AMP "%5s"   AMP "%5s"   AMP "%5s"   AMP "%6s"   AMP "%6s\n";
  char format2[256] ="xx2*| " "%6.4f" AMP "%7.5f" AMP "%5.2f" AMP "%5.3f" AMP "%5.3f" AMP "%5.3f" AMP "%5.3f" AMP "%5.3f" AMP "%5.3f" AMP "%6.4f" AMP "%6.4f\n";

  char format4[256] ="gg2*|  %5.3f&%6.4f &%5.1f & %5.2f   &%5.1f&%5.1f&%5.1f&%5.1f& %5.1f & %5.1f & %5.1f & %5.1f & %5.1f & %5.1f & %5.1f \n";
  char format5[256] ="mm1*| " "%6.4f"  " "  "%6.4f"  " "  "%6.4f"  " "  "%7.5f"  " "  "%5.2f"  " "  "%5.3f"  " "  "%5.3f"  " "  "%5.3f"  " "  "%5.3f"  " "  "%5.3f"  " "  "%5.3f"  " "  "%5.3f"  " "  "%5.3f"  " "  "%6.4f"  " "  "%6.4f" "\n" ;
  char format6[256] ="mmA*| %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n";


  /* BEGIN RNSID */

  /* EQUILIBRIUM VARIABLES */
  
  
  int    n_tab,                   /* Number of points in EOS file */
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
    **rho,                          /* potential \rho */ 
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

  /* =============================================== */
  /* Intermediate quantities used for output purpose */ 
  /* =============================================== */
  double Pa,Pe,td,beta,M_over_R;
  double invPa,invPe,invtd;
  /* FILE *file_2D;  */

 /* ============================================= */
 /*  Variable to save in TXT file                 */ 
 /* ============================================= */
 int varlistLEN = 15;
 struct {char name[128];double * valP;} varlist [] = {
      { "rhoc"    , &rho0_center} ,
      { "ec"      , &e_center   } ,
      { "rp_re"   , &r_ratio    } ,
      { "A_diff"  , &A_diff     } ,
      { "Re"      , &R_e	      } ,
      { "M0"      , &Mass_0     } ,
      { "M"       , &Mass	      } ,
      { "W"       , &W	      } ,
      { "T"       , &T	      } ,
      { "J"       , &J	      } ,
      { "beta"    , &beta	      } ,
      { "M_R"     , &M_over_R   } ,
      { "invPa"   , &invPa      } ,
      { "invPe"   , &invPe      } ,
      { "invtd"   , &invtd      } 
 };
  

  /* ******************************************
     .... a la charte
     ********************************************* */
  int Verbose = 0;
  int opterr,c;
  opterr = 0;
  for(c = getopt(argc,argv,"+hvur:a:t:f:k:g:d:e:x:y:z:l:m:n:i:o:s:");
      c != EOF;
      c = getopt(argc,argv,"+hvur:a:t:f:k:g:d:e:x:y:z:l:m:n:i:o:s:")) {
    switch(c) {
      /* ================================== */
      /* ===  CENTRAL DENSITY          ==== */
      /* ================================== */
    case 'r': /* Central density */
      sscanf(optarg,"%lf",&rho0_center);
      max_rho0_center   = rho0_center;
      min_rho0_center   = rho0_center;
      delta_rho0_center = 0.1e-3;
      break;
      /* SCANNING DIFFERENT VALUES OF THE CENTRAL DENSITY */
    case 'l':
      sscanf(optarg,"%lf",&min_rho0_center);
      break;
    case 'm':
      sscanf(optarg,"%lf",&max_rho0_center);
      break;
    case 'n':
      sscanf(optarg,"%lf",&delta_rho0_center);
      break;
      /* ================================== */
      /* ===  AXIS RATIO               ==== */
      /* ================================== */
    case 'a': /* axes ratio r_p/r_e */
      sscanf(optarg,"%lf",&r_ratio);
	  max_a_ratio = r_ratio;
	  min_a_ratio = r_ratio;
      break;
      /* SCANNING DIFFERENT VALUES OF THE AXIS RATIO */
    case 'x': /* axes ratio r_p/r_e MAX VALUE*/
      sscanf(optarg,"%lf",&max_a_ratio);
      break;
    case 'y': /* axes ratio r_p/r_e MIN_VALUE*/
      sscanf(optarg,"%lf",&min_a_ratio);
      break;
    case 'z': /* axes ratio r_p/r_e step*/
      sscanf(optarg,"%lf",&delta_a_ratio);
      break;
      /* ================================== */
      /* ===  EOS selection            ==== */
      /* ================================== */
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
      /* ================================== */
      /* ===  KIND OF ROTATION         ==== */
      /* ================================== */
    case 'd': /* \hat{A} diff rot law par */
      sscanf(optarg,"%lf",&A_diff);
	  printf("\n RNS: DIFFERENTIAL rotation sequence ! A = %g\n\n",A_diff);
      break;
    case 'u':
      strcpy(rotation_type,"uniform");
	  printf("\n RNS: Uniform rotation sequence !\n\n");
      break;
      /* ================================== */
      /* ===  OTHER SETTINGS           ==== */
      /* ================================== */
    case 'e': /* Accuracy goal */
      sscanf(optarg,"%lf",&accuracy);
      break;
    case 'i' : /* initial_model_data.h5  */
        sscanf(optarg,"%s",model2D_file_IN);
        if (strlen(model2D_file_IN) > 0 ) {
	  strcpy(recover_2Dmodel,"yes");
        }
        break;
    case 'o' : /* final_model_data.h5 */
      sscanf(optarg,"%s",model2D_file_OUT);
      break;
    case 's' : /* final_model_data.h5 */
      sscanf(optarg,"%s",export_data_file_name);
      break;
    case 'v':
      Verbose = 1;
      break;
    case 'h':
      printf("\n");
      printf("Help: ./RNS -options<meaning> {default value} \n");
      printf("\n");
      printf("      -r <central density>     {1.28e-3}\n");
      printf("      -a <axes ratio rp/re>    {1.0}\n");
      printf("      -x max(rp/re) -y nin(rp/re) -z delta(rp/re) \n");
      printf("      -l min(rhoc)  -m max(rhoc)  -n delta(rhoc) \n");
      printf("\n");
      printf("      -t <EOS type: poly/tab>  {poly} \n");
      printf("      -f <tab EOS file>        {no file} \n");
      printf("      -k <poly EOS K>          {100.0}\n");
      printf("      -g <poly EOS Gamma>      {2.0}\n");
      printf("\n");
      printf("      -d <{A} diff>            {1.0}\n");
      printf("      -u  (set uniform rotation)\n");
      printf("\n");
      printf("      -e<accuracy goal>       {1e-7}\n");
      printf("\n");
      printf("      -i <initial_model_data.hdf> \n");
      printf("      -o <final_model_data.hdf> \n");
      printf("      -s <export_data_file.txt> \n");
      printf("\n");
      return 0;
    default:
      break;
    }
  }
  
  /* ******************************************
     call rnsid
     ********************************************* */

  //if (Verbose == 1 ) {
  //    printf("RNSID called with the following parameters:  Verbose = True\n");
  //    printf("             RHO   = %16g %16g %16g \n",min_rho0_center,max_rho0_center,delta_rho0_center);
  //    printf("             rp_re = %16g %16g %16g \n",min_a_ratio,max_a_ratio,delta_a_ratio);
  //    printf("             A_diff= %16g \n",A_diff);
  //}
  
  /* Read valus from tabulated EOS and convert to C.U. */  

  if(strcmp(eos_type,"tab")==0) {
    printf(" TAB eos from file: %s\n",eos_file);    
    load_eos( eos_file, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, &n_tab );
  }  
  /* SET UP GRID */
  make_grid(s_gp, mu);
  /* ALLLOCATE MEMORY */
  rho = array_allocate(1,SDIV,1,MDIV);
  gama = array_allocate(1,SDIV,1,MDIV);
  alpha = array_allocate(1,SDIV,1,MDIV);
  omega = array_allocate(1,SDIV,1,MDIV);
  energy = array_allocate(1,SDIV,1,MDIV);
  pressure = array_allocate(1,SDIV,1,MDIV);
  enthalpy = array_allocate(1,SDIV,1,MDIV);
  velocity_sq = array_allocate(1,SDIV,1,MDIV);
  Omega_diff = array_allocate(1,SDIV,1,MDIV);


    /*==================================================
        EXPORT DATA IN TXT FILE
      ================================================== */ 
    if ( strlen(export_data_file_name) > 0 ) {
      export_data_file = fopen(export_data_file_name,"w");
      fprintf(export_data_file,"#######################################################\n");
      fprintf(export_data_file,"# Data file created by RNS\n");
      fprintf(export_data_file,"# ");
      for (ii=0; ii< argc ; ii++) 
        fprintf(export_data_file,"%s ",argv[ii]);
      fprintf(export_data_file,"\n");
      fprintf(export_data_file,"#\n");
      fprintf(export_data_file,"# column format: ");
      for (ii=0; ii< varlistLEN ; ii++) 
	fprintf(export_data_file, "%d:%s ",ii+1,(varlist[ii].name));
      fprintf(export_data_file,"\n#######################################################\n");
      fclose(export_data_file);
    }


 
  /* ============================================================= */
  /*       SCAN FOR DIFFERENT VALUES OF CENTRAL RHO                */
  /* ============================================================= */
  for(rho0_centerITR =min_rho0_center; rho0_centerITR <=  max_rho0_center; rho0_centerITR += delta_rho0_center) {
  rho0_center = rho0_centerITR;


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
  printf("e_center = %g\n",e_center);
 

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
  /* MAKE e_center DIMENSIONLESS FOR TAB. EOS */
  //
  // if(strcmp(eos_type,"tab")==0) 
  //   e_center *= (C*C*KSCALE);

  Omega_e=0.0; /* initialize ang. vel. at equator for diff. rot. */
  
  print_dif=1;
  
  /* COMPUTE DIMENSIONLESS CENTRAL PRESSURE AND ENTHALPY */

  make_center( e_center, log_e_tab, log_p_tab, log_h_tab, n_tab, 
	       eos_type, eos_k, Gamma_P, &p_center, &h_center); 
  
  /*what is this?*/
  //if(strcmp(eos_type,"tab")==0) 
  //      e_center *= (C*C*KSCALE);

  printf("Value of rho center is %16.8e\n",rho0_center);
  rho0_center =  (e_center+p_center)*exp(-h_center);
  printf("Value of rho center is %16.8e\n",rho0_center);
  
  if(print_dif==1){
    printf(" ****************************************************\n");
    printf(" ****************************************************\n");
    printf(" **                     RNSID                      **\n");
    printf(" **      ROTATING NEUTRON STAR INITIAL DATA        **\n");
    printf(" **                                                **\n");
    printf(" **      SEQUENCE OF DIFFERENT AXES RATIO          **\n");
    printf(" **                                                **\n");
    printf(" ****************************************************\n");
    printf(" ****************************************************\n");
  }   
  

  if ( strcmp(recover_2Dmodel, "no") == 0 ) {          
     printf(" Iterating equilibrium model \n");
 
     /* COMPUTE A SPHERICAL STAR AS A FIRST GUESS */
     printf(" Guessing with sphere (RNSID) \n");
     guess( s_gp, eos_type, eos_k, e_center, p_center, p_surface, e_surface, 
	     Gamma_P, log_e_tab, log_p_tab, log_h_tab, n_tab, rho, gama, 
	     alpha, omega, &r_e );     


    /* ========================================================================
       First cicle to reduce r_ratio close to the initial one using step -0.1  
       ========================================================================= */
    printf(" Lowering a_ratio\n");
    for(r_ratio=1.0;r_ratio > max_a_ratio; r_ratio -=0.1) {
    	printf(" Reducing a_ratio: Iteration for a_ratio = %g \n",r_ratio);
    	iterate( s_gp, mu, eos_type, eos_k, log_e_tab, log_p_tab, log_h_tab, 
    	         n_tab, Gamma_P, r_ratio, h_center, enthalpy_min, a_check, 
    	         accuracy, print_dif, cf, &r_e, rho, gama, alpha, omega, 
    	         energy, pressure, enthalpy, velocity_sq, &Omega,
    	         rotation_type, A_diff, &Omega_e, Omega_diff, RNS_lmax);
    }
  } else {
    hdf5_read_var(&sdiv,&mdiv,model2D_file_IN,eos_type,eos_file,&eos_k,&Gamma_P,rotation_type,&A_diff,&r_ratio,&rho0_center,&r_e,
                  s_gp, mu,rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq, &Omega,&Omega_e,Omega_diff);
    comp_values( s_gp, mu, r_ratio, e_surface, r_e, eos_type, log_e_tab,
    		   log_n0_tab, n_tab, Omega, rho, gama, alpha, omega, 
    	       energy, pressure, enthalpy, velocity_sq, &Mass, 
    	       &Mass_0, &T, &W, &Omega_K, &R_e, rotation_type, Omega_diff,
    	       &J);

  }

  /* ========================================================================
     Main cicle over r_ratio using step -(delta_a_ratio)  
  ========================================================================= */
  for(r_ratio=max_a_ratio;r_ratio+0.001 >min_a_ratio; r_ratio -=delta_a_ratio) {
    idx +=1;
    printf(" Computing (%d) iteration for rho_c = %g , a_ratio = %g \n",idx,rho0_center,r_ratio);
  
    iterate( s_gp, mu, eos_type, eos_k, log_e_tab, log_p_tab, log_h_tab, 
             n_tab, Gamma_P, r_ratio, h_center, enthalpy_min, a_check, 
             accuracy, print_dif, cf, &r_e, rho, gama, alpha, omega, 
             energy, pressure, enthalpy, velocity_sq, &Omega,
             rotation_type, A_diff, &Omega_e, Omega_diff, RNS_lmax);
    
    /* COMPUTE EQUILIBRIUM QUANTITIES (Mass, Radius, T/W etc.) */
    comp_values( s_gp, mu, r_ratio, e_surface, r_e, eos_type, log_e_tab,
		   log_n0_tab, n_tab, Omega, rho, gama, alpha, omega, 
	       energy, pressure, enthalpy, velocity_sq, &Mass, 
	       &Mass_0, &T, &W, &Omega_K, &R_e, rotation_type, Omega_diff,
	       &J);
  
  
    /* RETURN OMEGA AND R_E VALUES */
    (*Omega_pt) = Omega;
	(*R_e_pt) = R_e;
    (*r_e_pt) = r_e;
	(*mass0_pt) = Mass_0;
  

    if (strcmp(rotation_type,"uniform")==0) {
        invPa = Omega                  * (cactusT/(1000*2*PI));
        invPe = Omega                  * (cactusT/(1000*2*PI));
    } else { 
        invPa = Omega_diff[1][1]      * (cactusT/(1000*2*PI));
        invPe = Omega_e               * (cactusT/(1000*2*PI));
    }
    invtd = sqrt(R_e*R_e*R_e/Mass_0) * (cactusT/(1000*2*PI));
    beta     = T/W;
    M_over_R = Mass/R_e;

   /* PRINT-OUT SOME EQUILIBRIUM QUANTITIES */
   if(print_dif==1){   
     printf("Equilibrium model done in c=G=M_sun=1 dimensionless form (%s)\n",rotation_type); 
     printf("===================================================================\n");
     printf(" rho_center   e_center    Mass      Mass_0      R_e\n");
     printf(" %5.4e %5.4e %5.4e %5.4e %5.4e\n", rho0_center, e_center, Mass, Mass_0, R_e);  
     if (strcmp(rotation_type,"uniform")==0) {
       printf("     J         T/W       Omega   Omega_Kepler axes_ratio  J/M^2 \n");
       printf(" %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e\n",  ( (Omega > 0.0) ? J :0.0), T/W,
     	  Omega, Omega_K, r_ratio,J/(Mass*Mass));
     } else {
       printf("     J         T/W       Omega_c Omega_e   Omega_Kepler axes_ratio  J/M^2 \n");
       printf(" %5.4e %5.4e %5.4e %5.4e %5.4e  %5.4e %5.4e\n",  ( (Omega > 0.0) ? J : 0.0), T/W,
     	  Omega, Omega_e, Omega_K,r_ratio,J/(Mass*Mass));
     }
     printf("===================================================================\n");
     
     if (Verbose == 1 ) {
         printf("yy| rhoc(%d)   = %16.8e ;\n",idx,rho0_center);
         printf("yy| ec(%d)     = %16.8e ;\n",idx,e_center);
         printf("yy| rp_re(%d)  = %16.8e ;\n",idx,r_ratio);
         printf("yy| A_diff(%d) = %16.8e ;\n",idx,A_diff);
         printf("yy| Re(%d)     = %16.8e ;\n",idx,R_e);
         printf("yy| M0(%d)     = %16.8e ;\n",idx,Mass_0);
         printf("yy| M(%d)      = %16.8e ;\n",idx,Mass);
         printf("yy| W(%d)      = %16.8e ;\n",idx,W);
         printf("yy| T(%d)      = %16.8e ;\n",idx,T);
         printf("yy| J(%d)      = %16.8e ;\n",idx,J);
         printf("yy| beta(%d)   = %16.8e ;\n",idx,beta);
         printf("yy| M_R(%d)    = %16.8e ;\n",idx,M_over_R);
         printf("  | === periods in msecs ==========\n");
         printf("yy| invPa(%d)   = %16.8e ;\n",idx,invPa);
         printf("yy| invPe(%d)   = %16.8e ;\n",idx,invPe);
         printf("yy| invtd(%d)   = %16.8e ;\n",idx,invtd);
         printf("yy| invPK(%d)   = %16.8e ;\n",idx,Omega_K  * (cactusT/(1000*2*PI)));
         printf("yy| MDIV = %d ; SDIV = %d;\n",MDIV,SDIV);
     
         printf("xx0*|============================================================================================\n");	
         printf(format1,"rho_c","r_p/r_e","R_e","M_0","M","J","P_a" ,"P_e" ,"t_D"  ,"M/R_e","beta");
         printf(format1,"(e-4)",""       ,""   ,""   ,"" ,"" ,"msec","msec","msec" ,""     ,"");
         printf("xx0*|--------------------------------------------------------------------------------------------\n");	
         printf(format2,10000.0*rho0_center,r_ratio,R_e,Mass_0,Mass,J,1.0/invPa,1.0/invPe,1.0/invtd,M_over_R,beta);
         printf("xx3*|============================================================================================\n");
    }

    /*==================================================
        EXPORT DATA IN TXT FILE
      ================================================== */ 
    if ( strlen(export_data_file_name) > 0 ) {
      /* printf("%s %16e ",export_data_file_name,1.0); */
      export_data_file = fopen(export_data_file_name,"a");
      for (ii=0; ii< varlistLEN ; ii++) 
	fprintf(export_data_file,"%.16e ",*(varlist[ii].valP));
      fprintf(export_data_file,"\n");
      fclose(export_data_file);
    }

	/* =====================================================================
    printf("mm0*|==============================================================================================================================\n");
    printf(format5,M_over_R,beta,10000.0*rho0_center,r_ratio,R_e,Mass_0,Mass,J,J/(Mass*Mass),Pa,Pe,Pe/Pa,td,100.0*T,100.0*W);
    printf(format6,M_over_R,beta, rho0_center,r_ratio,R_e,Mass_0,Mass,J,J/(Mass*Mass),Pa,Pe,Pe/Pa,td,T,W);
    printf("mm2*|==============================================================================================================================\n");

    printf("gg0*|=================================================================================================================\n");
    printf("gg1*| M/R_e & beta  &rho_c & r_p/r_e & R_e & M_0 & M   & J   & J/M^2 & P_a   & P_e   & Pe/Pa & t_D   &  T    &  W     \n");
    printf("gg1*|       &       &(e-4) &         &     &     &     &     &       & msec  & msec  &       & msec  & (e-2) & (e-2)  \n");
    printf("gg0*|-----------------------------------------------------------------------------------------------------------------\n");
    printf(format4,M_over_R,beta,10000.0*rho0_center,r_ratio,R_e,Mass_0,Mass,J,J/(Mass*Mass),Pa,Pe,Pe/Pa,td,100.0*T,100.0*W);
    printf("gg3*|=================================================================================================================\n");

	====================================================================== */

   }


  }
  r_ratio +=delta_a_ratio;
  /* ========================================================================
     AND OF Main cicle over r_ratio using step -(delta_a_ratio)  
     for(r_ratio=max_a_ratio;r_ratio+0.001 >min_a_ratio; r_ratio -=delta_a_ratio) 
     CORRECTING r_ratio to the last used in computation
  ========================================================================= */

  }
  /* ============================================================= */
  /*       END OF SCAN FOR DIFFERENT VALUES OF CENTRAL RHO         */
  /* ============================================================= */


  /* ===========================================================================
   Final step to save the OUTPUT of the model as HDF5 File (if required) 
   ========================================================================== */ 
  int lenmodel2D_file_OUT=strlen(model2D_file_OUT);
  int sdiv1,mdiv1;
  char eos_type1[128],eos_file1[128],rotation_type1[128];
  double A_diff1,r_ratio1,rho0_center1;
  double eos_k1,Gamma_P1;
  double **rho1;

  if (lenmodel2D_file_OUT > 0) {
    printf("Length of model2D_file_OUT is %d (%s)\n",lenmodel2D_file_OUT,model2D_file_OUT);
    hdf5_save_var(&sdiv,&mdiv,model2D_file_OUT,
                  eos_type,eos_file,
                  &eos_k,&Gamma_P,
                  rotation_type,&A_diff,
                  &r_ratio,&rho0_center,&r_e,
                  s_gp, mu,
                  rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq, &Omega,&Omega_e,Omega_diff);
  
    rho1 = array_allocate(1,SDIV,1,MDIV);
  
    hdf5_read_var(&sdiv1,&mdiv1,model2D_file_OUT,eos_type1,eos_file1,&eos_k1,&Gamma_P1,rotation_type1,&A_diff1,&r_ratio1,&rho0_center1,&r_e,
                  s_gp, mu,rho1, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq, &Omega,&Omega_e,Omega_diff);
    printf("%d %d   --- %d %d \n",sdiv,mdiv,sdiv1,mdiv1);
    printf("%f %f %f   --- %f %f %f \n",A_diff,r_ratio,rho0_center,A_diff1,r_ratio1,rho0_center1);
    printf("%f %f     --- %f %f \n",eos_k,Gamma_P,eos_k1,Gamma_P1);
    printf("(%s)-(%s)-(%s) --- (%s)-(%s)-(%s) \n",rotation_type,eos_type,eos_file,rotation_type1,eos_type1,eos_file1);
    for (i = 0; i < SDIV; i++) {
      for (j = 0; j < MDIV; j++) {
    	if ( rho[i+1][j+1] != rho1[i+1][j+1] ) {
    	  printf("elements %d,%d of rho are different %f != %f\n",i,j,rho[i+1][j+1], rho1[i+1][j+1]);
        }
      }
    }
    //int sdiv=SDIV;
    //int mdiv=MDIV;
    //hdf5_read_var(&sdiv,&mdiv,model2D_file_OUT,rho, gama, alpha, omega, 
    //		  energy, pressure, enthalpy, velocity_sq, &Omega,
    //              rotation_type, A_diff, &Omega_e,Omega_diff);
  
    hdf5_read_var(&sdiv,&mdiv,model2D_file_OUT,eos_type,eos_file,&eos_k,&Gamma_P,rotation_type,&A_diff,&r_ratio,&rho0_center,&r_e,
                  s_gp, mu,rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq, &Omega,&Omega_e,Omega_diff);
    comp_values( s_gp, mu, r_ratio, e_surface, r_e, eos_type, log_e_tab,
    		   log_n0_tab, n_tab, Omega, rho, gama, alpha, omega, 
    	       energy, pressure, enthalpy, velocity_sq, &Mass, 
    	       &Mass_0, &T, &W, &Omega_K, &R_e, rotation_type, Omega_diff,
    	       &J);
    invPa = Omega_diff[1][1]         * (cactusT/(1000*2*PI));
    invPe = Omega_e                  * (cactusT/(1000*2*PI));
    invtd = sqrt(R_e*R_e*R_e/Mass_0) * (cactusT/(1000*2*PI));
    beta     = T/W;
    M_over_R = Mass/R_e;
    printf("xx0*|============================================================================================\n");	
    printf("xx0*|             Properties of the save model                                                   \n");	
    printf("xx0*|============================================================================================\n");	
    printf(format1,"rho_c","r_p/r_e","R_e","M_0","M","J","P_a" ,"P_e" ,"t_D"  ,"M/R_e","beta");
    printf(format1,"(e-4)",""       ,""   ,""   ,"" ,"" ,"msec","msec","msec" ,""     ,"");
    printf("xx0*|--------------------------------------------------------------------------------------------\n");	
    printf(format2,10000.0*rho0_center,r_ratio,R_e,Mass_0,Mass,J,1.0/invPa,1.0/invPe,1.0/invtd,M_over_R,beta);
    printf("xx3*|============================================================================================\n");
  }
  
  
  return 0;
}
