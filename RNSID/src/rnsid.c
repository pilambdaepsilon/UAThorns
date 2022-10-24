/**************************************************************************
*                              RNSID.C                                    *
*                                                                         *
*                             FOR CACTUS                                  *
*                                                                         *
*  Author:  Nikolaos Stergioulas                                          *
*                                                                         *
*  E-mail: niksterg@astro.auth.gr                                         *
*                                                                         *
*                                                                         *
*  First Release Date:    October 1998                                    *
*                                                                         *
**************************************************************************/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "consts.h"
#include "equil.h"
#include "equil_util.h"
#include "rnsid_util.h"
#include "rnsid.h"
#include "hdf5_save.h"


#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sx (&scon[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sy (&scon[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sz (&scon[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

static inline int exists_file_name(const char *fname)
{
	FILE *file;
	if (file = fopen(fname, "r"))
	{
	    fclose(file);
	    return 1;
	}
	return 0;
}

/*************************************************************************/
/* COMPUTE CARTESIAN INITIAL DATA FOR A ROTATING NEUTRON STAR            */
/*************************************************************************/
void hydro_rnsid(const cGH *cctkGH,
          CCTK_REAL *x_grid,
          CCTK_REAL *y_grid,
          CCTK_REAL *z_grid,
          CCTK_REAL eos_k,
          CCTK_REAL eos_ideal_fluid_gamma,
          CCTK_REAL rnsid_rho_min,
          CCTK_REAL *Omega_pt,
          CCTK_REAL *R_e_pt,
          CCTK_REAL *r_e_pt,
          CCTK_REAL *mass0_pt,
          CCTK_REAL *gamma_center)


{ /* BEGIN RNSID */

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

 /* EQUILIBRIUM VARIABLES */

 int    n_tab,                        /* Number of points in EOS file */
        n_tab_beta,
        a_check=0,                    /* if =200, iteration diverges */
        print_dif,                    /* if =1, monitor convergence */
        i;                            /* loop counter */


 double log_e_tab[MAX_NTAB],          /* energy dens./c^2 in tab. EOS */
        log_p_tab[MAX_NTAB],          /* pressure in tabulated EOS */
        log_h_tab[MAX_NTAB],          /* enthalpy in EOS file */
        log_n0_tab[MAX_NTAB],         /* number density in EOS file */
        log_e_tab_beta[MAX_NTAB],
        Y_e_tab[MAX_NTAB],
        Gamma_tab[MAX_NTAB],          /* Gamma in tab. EOS file */
        e_center,                     /* central energy density */
        Y_e_center,
        p_center,                     /* central pressure */
        h_center,                     /* central enthalpy */
        e_surface,                    /* surface en. density */
        p_surface,                    /* surface pressure */
        enthalpy_min,                 /* minimum enthalpy in EOS */
        *s_gp,                        /* s grid points */
        *mu,                          /* \mu grid points */
      **rho_potential,                          /* potential \rho_potential */
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

 //  char eos_file[256] = "no EOS file specified",   /* EOS file name */
 //        eos_type[80]  = "poly";                    /* EOS type (poly or tab) */
 char * message;

 /* INITIAL DATA VARIABLES */

 int m,                               /* counter */
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
       Y_e_ijk,
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
       h_ijk,                         /* h = 1 + eps + P/rho_potential */
       distance_ijk =0,               /* Signed distance to surface */
       e_atm,                         /* energy density of atmosphere */
       p_atm,                         /* pressure of atmosphere */
       Y_e_atm,
       e_atm_beta,
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

 FILE *file_2D;

 int nx=cctkGH->cctk_lsh[0]; int ny=cctkGH->cctk_lsh[1]; int nz=cctkGH->cctk_lsh[2];
int n_nearest_beta;

 CCTK_REAL rho0_center;

 /*

 HISTORICAL NOTE ON NAMES OF VARIABLES:

 old name               new name (from version 1.25 on)

 pert_amp               pert_amplitude (now a parameter)
 Gamma_P                eos_ideal_fluid_gamma (now a parameter)
 r_ratio                axes_ratio (now a parameter)
 rho                    rho_potential

 */


      /* COMPUTE POLYTROPIC INDEX AND CENTRAL ENERGY DENSITY */
      n_P=1.0/(eos_ideal_fluid_gamma-1.0);
      e_center = (eos_k*pow(rho_central,eos_ideal_fluid_gamma)/(eos_ideal_fluid_gamma-1.0)+rho_central);

      /* TABULATED EOS OPTION */
      if(strcmp(eos_type,"tab")==0) {
        /* --V0-- load_eos( eos_file, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, Gamma_tab, &n_tab ); */
        int n_nearest;
        double n0;
        /* ==================================================== */
        /* printf(" TAB eos from file: %s\n",eos_file); */
        message = (char *)malloc(200*sizeof(char));
        sprintf(message," TAB eos from file: %s",eos_file);
        CCTK_INFO(message);
        free(message);
        /* ==================================================== */
        load_eos( eos_file, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, &n_tab );
	n_nearest = 50;
	n0 = rho_central/(MB*cactusM);
	e_center = pow(10.0,interp(log_n0_tab, log_e_tab, n_tab,log10(n0), &n_nearest));
      }

      if(use_beta_equil_Y_e){
        message = (char *)malloc(200*sizeof(char));
        sprintf(message, "Reading beta_equilibrium file for Y_e initialization %s",beta_equil_file);
        CCTK_INFO(message);
        free(message);
        /* ==================================================== */
        load_beta_equil( beta_equil_file, log_e_tab_beta, Y_e_tab, &n_tab_beta );
        Y_e_center = interp(log_e_tab_beta, Y_e_tab, n_tab_beta,log10(e_center), &n_nearest_beta);
        Y_e_atm = Y_e_tab[1];
        e_atm_beta = pow(10.0,log_e_tab_beta[1]);
        CCTK_VInfo(CCTK_THORNSTRING, "WUTANG:  [central Y_e]%1.5e   [ntab]%d   [e_center-Cactus]%1.5e    [e_center]%1.5e    [lowest e]%1.5e    [Y_e_atm]%1.5e", Y_e_center, n_tab_beta, e_center, e_center/(cactusM*cactusV), e_atm_beta/(cactusM*cactusV), Y_e_atm);
      } 


      /* SET UP GRID */
      s_gp=malloc((SDIV+1)*sizeof(double));
      mu=malloc((MDIV+1)*sizeof(double));
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

      /* INITIALIZE VARIABLES WITH ZERO */

      #pragma omp parallel for
      for(s=1;s<=SDIV;s++)
        for(m=1;m<=MDIV;m++) {
          rho_potential[s][m] = 0.0e0;
          gama[s][m] = 0.0e0;
          alpha[s][m] = 0.0e0;
          omega[s][m] = 0.0e0;
          energy[s][m] = 0.0e0;
          pressure[s][m] = 0.0e0;          
	        enthalpy[s][m] = 0.0e0;
          velocity_sq[s][m] = 0.0e0;
          Omega_diff[s][m] = 0.0e0;
        }



      /* SET DEFAULT EQUILIBRIUM PARAMETERS */

      if(strcmp(eos_type,"tab")==0) {
        e_surface=7.8e-15;
        p_surface=1.12379e-28;
        enthalpy_min=1.0/(C*C);
      }
      else {
            e_surface=0.0;
            p_surface=0.0;
            enthalpy_min=0.0;
      }

      Omega_e=0.0; /* initialize ang. vel. at equator for diff. rot. */

      print_dif=1;


      z_print = nz/4;

 
      /* CHANGE e_center TO POLYTROPIC DIMENSIONLESS UNITS */
      /*      e_center /= ( 1.0/pow(eos_k, n_P) ); */

      /* MAKE e_center DIMENSIONLESS FOR TAB. EOS */
      //
      // if(strcmp(eos_type,"tab")==0) 
      //   e_center *= (C*C*KSCALE);


      /* COMPUTE DIMENSIONLESS CENTRAL PRESSURE AND ENTHALPY */
 
      /*-V0-- make_center( e_center, log_e_tab, log_p_tab, log_h_tab, n_tab,      */
      /*-V0--             eos_type, eos_ideal_fluid_gamma, &p_center, &h_center); */
      make_center( e_center, log_e_tab, log_p_tab, log_h_tab, n_tab, 
                   eos_type, eos_k,eos_ideal_fluid_gamma, &p_center, &h_center); 

      rho0_center =  (e_center+p_center)*exp(-h_center);



      if(print_dif==1){
        /* ====================================================================
        printf(" ****************************************************\n");
        printf(" ****************************************************\n");
        printf(" **              HYDRO - RNSID                     **\n");
        printf(" **      ROTATING NEUTRON STAR INITIAL DATA        **\n");
        printf(" ****************************************************\n");
        printf(" ****************************************************\n");
        if( strcmp(recover_2Dmodel, "no")==0)         
	  printf(" Iterating equilibrium model\n");
	  ======================================================================*/
        CCTK_INFO(" ****************************************************");
        CCTK_INFO(" ****************************************************");
        CCTK_INFO(" **                   RNSID                        **");
        CCTK_INFO(" **      ROTATING NEUTRON STAR INITIAL DATA        **");
        CCTK_INFO(" ****************************************************");
        CCTK_INFO(" ****************************************************");
        if( strcmp(recover_2Dmodel, "no")==0)         
	  CCTK_INFO(" Iterating equilibrium model");
      }   


      /* EITHER COMPUTE THE MODEL AND PROCEED OR COMPUTE IT AND SAVE IT
	 IN A FILE */
      
      int recovered_from_2d = 0;
      if(  strcmp(recover_2Dmodel, "yes")==0 && exists_file_name(model2D_file)==1) {

             /* RECOVER FROM 2D FILE */
             /* ==================================================== */
             /* printf(" Recovering 2D model form file %s\n",model2D_file); */
             message = (char *)malloc(200*sizeof(char));
             sprintf(message," Recovering 2D model form file %s",model2D_file);
             CCTK_INFO(message);
             free(message);
             /* ==================================================== */
             
	     int sdiv,mdiv;
             hdf5_read_var(&sdiv,&mdiv,model2D_file,
                           eos_type,eos_file,&eos_k,&eos_ideal_fluid_gamma,rotation_type,
                           &A_diff,&axes_ratio,&rho0_center,&r_e,
                           s_gp,mu,rho_potential, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq, &Omega,&Omega_e,Omega_diff);
             recovered_from_2d = 1;
      } /* END RECOVER FROM 2D FILE, so also end of finding the equil state */
      else {
 
	  printf(" Guessing with a sphere (rnsid.c) [p_surface]: %e", p_surface);
          ///* COMPUTE A SPHERICAL STAR AS A FIRST GUESS */
          guess( s_gp, eos_type, eos_k,e_center, p_center, p_surface, e_surface, 
	       eos_ideal_fluid_gamma, log_e_tab, log_p_tab, log_h_tab, n_tab, rho_potential, gama, 
               alpha, omega, &r_e );     
	
          /* If the axis ratio is less than 0.8 (0.6), one needs to compute 
             the model with axes_ratio=0.8 (0.6) first and then compute the 
             desired model (otherwise the iteration may not converge). The 
             subroutine "iterate" starts with the current guess and 
             converges to the desired rotating model. The guess is either a 
             spherical star or a previously computed (slower) rotating 
             star. 
          */
          /* ========================================================================
              First cicle to reduce r_ratio close to the initial one using step -0.1  
            ========================================================================= */
          double r_ratio;
          for(r_ratio=1.0;r_ratio > axes_ratio; r_ratio -=0.1) {
	      //printf(" Reducing a_ratio: Iteration for a_ratio = %g \n",r_ratio);
              message = (char *)malloc(200*sizeof(char));
	      sprintf(message," Reducing a_ratio: Iteration for a_ratio = %g",r_ratio);
              CCTK_INFO(message);
              free(message);
     	      iterate( s_gp, mu, eos_type, eos_k, log_e_tab, log_p_tab, log_h_tab, 
    	               n_tab, eos_ideal_fluid_gamma, r_ratio, h_center, enthalpy_min, a_check, 
    	               accuracy, print_dif, cf, &r_e, rho_potential, gama, alpha, omega, 
    	               energy, pressure, enthalpy, velocity_sq, &Omega,
    	               rotation_type, A_diff, &Omega_e, Omega_diff, RNS_lmax);
          }


	  // CCTK_REAL cf = 1.0;
	  CCTK_INFO(" calling iterate()");
          iterate( s_gp, mu, eos_type, eos_k, log_e_tab, log_p_tab, log_h_tab, 
                   n_tab, eos_ideal_fluid_gamma, axes_ratio, h_center, enthalpy_min, a_check, 
                   accuracy,print_dif,cf, &r_e, rho_potential, gama, alpha, omega, 
                   energy, pressure, enthalpy, velocity_sq, &Omega,
                   rotation_type,A_diff,&Omega_e, Omega_diff,RNS_lmax);

      } 

      if( strcmp(save_2Dmodel, "yes")==0 && recovered_from_2d == 0 ) {

        /* SAVE 2D FILE */
        int sdiv,mdiv;
        sdiv=SDIV;
        mdiv=MDIV;
        /* Only the root process saves */
        if (CCTK_MyProc(cctkGH) == 0) {
          hdf5_save_var(&sdiv,&mdiv,model2D_file,
                        eos_type,eos_file,&eos_k,&eos_ideal_fluid_gamma,rotation_type,
                        &A_diff,&axes_ratio,&rho0_center,&r_e,
                        s_gp,mu,rho_potential, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq, &Omega,&Omega_e,Omega_diff);
        }
	if (exists_file_name(model2D_file) == 0) {
	  message = (char *)malloc(200*sizeof(char));
	  sprintf(message,"File %s could not be written!",model2D_file);
	  CCTK_WARN(CCTK_WARN_PICKY, message);
	  free(message);
        }
      }  /* END SAVE 2D FILE */


      /* COMPUTE EQUILIBRIUM QUANTITIES (Mass, Radius, T/W etc.) */ 

      comp_values( s_gp, mu, axes_ratio, e_surface, r_e, eos_type, log_e_tab,
                   log_n0_tab, n_tab, Omega, rho_potential, gama, alpha, omega, 
                   energy, pressure, enthalpy, velocity_sq, &Mass, 
                   &Mass_0, &T, &W, &Omega_K, &R_e, rotation_type,Omega_diff,
                   &J);


      /* TRANSFORM UNITS TO c=G=M_sun=1 */
      /*
      transform_units( eos_type, n_P, eos_k, &rho0_center, &e_center, 
                       &p_center, &r_e, omega, energy, pressure, &Mass, 
                       &Mass_0, &T, &W, &Omega, &Omega_K, &R_e,
                       &Omega_e, Omega_diff, &J);
      */

      /* RETURN OMEGA AND R_E VALUES */

      (*Omega_pt) = Omega;
      (*R_e_pt) = R_e;
      (*r_e_pt) = r_e;
      (*mass0_pt) = Mass_0;


      /* PRINT-OUT SOME EQUILIBRIUM QUANTITIES */
      if (print_dif==1) {
        message = (char *)malloc(200*sizeof(char));
        sprintf(message," %5.4e %5.4e %5.4e %5.4e %5.4e", rho0_center, e_center, Mass, Mass_0, R_e);  
	CCTK_INFO("Equilibrium model done in c=G=M_sun=1 dimensionless form"); 
        CCTK_INFO(" rho_center   e_center    Mass      Mass_0      R_e");
        CCTK_INFO(message);
        if (strcmp(rotation_type,"uniform")==0) {
          CCTK_INFO("     J         T/W       Omega   Omega_Kepler axes_ratio  J/M^2 ");
          sprintf(message," %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e",  ( (Omega > 0.0) ? J :0.0), T/W,
		  Omega, Omega_K, axes_ratio,J/(Mass*Mass));
	  CCTK_INFO(message);  
        } else {
          CCTK_INFO("     J         T/W       Omega_c Omega_e   Omega_Kepler axes_ratio  J/M^2 ");
          sprintf(message," %5.4e %5.4e %5.4e %5.4e %5.4e  %5.4e %5.4e",  ( (Omega > 0.0) ? J : 0.0), T/W,
		  Omega, Omega_e, Omega_K,axes_ratio,J/(Mass*Mass));
	  CCTK_INFO(message);  
        }
        free(message);
      }

      /* ==============================================================================
      if(print_dif==1){
        printf("In c=G=M_sun=1 dimensionless form\n");
        printf(" Equilibrium model done\n");
        printf(" rho_center   e_center    Mass      Mass_0      R_e\n");  
        fflush(stdout);
        printf(" %5.4e %5.4e %5.4e %5.4e %5.4e\n", rho0_center, e_center,
                 Mass, Mass_0, R_e);
        fflush(stdout);
        if(strcmp(rotation_type,"uniform")==0) {
         printf("     J         T/W       Omega   Omega_Kepler axes_ratio\n");
         fflush(stdout);
         printf(" %5.4e %5.4e %5.4e %5.4e %5.4e\n",  ( (Omega > 0.0) ?
                                                       J :
                                                       0.0), T/W,
                 Omega, Omega_K, axes_ratio);  
         fflush(stdout);
         printf("    J/M^2\n");
         fflush(stdout);
         printf(" %5.4e\n",J/(Mass*Mass));
         fflush(stdout);
        }
        else {
          printf("     J         T/W       Omega_c    Omega_e    Omega_K\n");
          fflush(stdout);
          printf(" %5.4e %5.4e %5.4e %5.4e %5.4e\n",  ( (Omega > 0.0) ?
                                                        J :
                                                        0.0), T/W,
                 Omega, Omega_e, Omega_K);
          fflush(stdout);
            printf("    J/M^2   axes_ratio\n");
             fflush(stdout);
            printf(" %5.4e %5.4e\n",J/(Mass*Mass), axes_ratio);
          fflush(stdout);
         } 
      }
     ==================================================================================== */

    /* CONSTRUCT ARRAYS WITH NEEDED POLAR QUANTITIES */

      CCTK_INFO("constructing array with needed polar quantities");
      nu = array_allocate(1,SDIV,1,MDIV);      
      B = array_allocate(1,SDIV,1,MDIV);       
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

      CCTK_INFO("computing ID");
      rho_0_atm = rnsid_rho_min; /* rename the constant for historical reasons */
      e_atm = rho_0_atm;
      if(strcmp(eos_type,"tab")==0) {
        int n_nearest = n_tab/2;
        p_atm = p_at_e(e_atm, log_p_tab, log_p_tab, n_tab, &n_nearest);
      } else { if(strcmp(eos_type,"poly")==0) {
          p_atm = eos_k*pow(rho_0_atm,eos_ideal_fluid_gamma);
        }
      }
      CCTK_INFO("starting loop");

#pragma omp for collapse (3)			\
  private(i,j,k) \
  private(nu_ijk,B_ijk,alpha_ijk,omega_ijk,nu_dr_ijk,B_dr_ijk, \
          nu_dtheta_ijk, B_dtheta_ijk,alpha_dtheta_ijk,    \
          omega_dr_ijk, omega_dtheta_ijk, \
          rho_0_ijk, energy_ijk, pressure_ijk,distance_ijk,Omega_ijk, \
          exp_nu_ijk,exp_alpha_ijk,r_bar_ijk, \
          dr_dx,dr_dy,dr_dz,dtheta_dx,dtheta_dy,dtheta_dz, \
          alpha_dx,alpha_dy,omega_dx,omega_dy,omega_dz) 
      for(i=1;i<=nx;i++)
         for(j=1;j<=ny;j++) 
           for(k=1;k<=nz;k++) {

               x_i = x_grid[i-1+nx*(j-1+ny*(k-1))];
               y_j = y_grid[i-1+nx*(j-1+ny*(k-1))];
               z_k = z_grid[i-1+nx*(j-1+ny*(k-1))];

               grid_interp_all( s_gp, mu, r_e, nx, ny, nz, 
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

	       if(use_beta_equil_Y_e){
                 Y_e_ijk = interp(log_e_tab_beta, Y_e_tab, n_tab_beta,log10(energy_ijk), &n_nearest_beta);
               }


              /* *************************************** */
              /* DETECT if it is in the ATMOSPHERE       */
              /* *************************************** */
              if( (rho_0_ijk<=0.0) || (energy_ijk<=0.0) || 
                  (pressure_ijk<=0.0) ) {
 
                rho_0_ijk   = rho_0_atm;
                energy_ijk  = e_atm + 1.e-20;
                pressure_ijk= p_atm;
		if(use_beta_equil_Y_e){Y_e_ijk = Y_e_atm;}

              }
	      if(use_beta_equil_Y_e){
                if(Y_e_ijk < 0.036){
                //  CCTK_VInfo(CCTK_THORNSTRING, "WUTANG YE SMALLER THAN LIMIT %1.3e %1.3e", Y_e_ijk, Y_e_atm);
                  Y_e_ijk = 0.036;
                }
              }


              /* *************************************** */
              /* END ATMOSPERE SETTINGs                      */
              /* *************************************** */

              exp_nu_ijk = exp(nu_ijk);
              exp_alpha_ijk = exp(alpha_ijk);

              r_ijk = sqrt(SQ(x_i)+SQ(y_j)+SQ(z_k));
              r_bar_ijk = sqrt(SQ(x_i)+SQ(y_j));


              alp[i-1+nx*(j-1+ny*(k-1))] = exp_nu_ijk;


              if(x_i==0.0 && y_j==0.0) {

                gxx[i-1+nx*(j-1+ny*(k-1))] = SQ(exp_alpha_ijk);
                gyy[i-1+nx*(j-1+ny*(k-1))] = SQ(exp_alpha_ijk);
                gzz[i-1+nx*(j-1+ny*(k-1))] = SQ(exp_alpha_ijk);

                gxy[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                gxz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                gyz[i-1+nx*(j-1+ny*(k-1))] = 0.0;

                kxx[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                kyy[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                kzz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                kxy[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                kxz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                kyz[i-1+nx*(j-1+ny*(k-1))] = 0.0;

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

                /* enforce omega_dz=0 at z=0 (it is slightly nonzero due
                   to O(h) forwards formula in computing derivative) */

                if(z_k==0.0)
                  omega_dz = 0.0;
                else
                  omega_dz = dr_dz*omega_dr_ijk + dtheta_dz*omega_dtheta_ijk;


                gxx[i-1+nx*(j-1+ny*(k-1))] = ( SQ(B_ijk*y_j/exp_nu_ijk)
                                             +SQ(exp_alpha_ijk*x_i) ) /
                                             (SQ(x_i)+SQ(y_j));

                gxy[i-1+nx*(j-1+ny*(k-1))] = ( SQ(exp_alpha_ijk)
                                            -SQ(B_ijk/exp_nu_ijk) ) *
                                            x_i*y_j/(SQ(x_i)+SQ(y_j));

                gxz[i-1+nx*(j-1+ny*(k-1))] = 0.0;

                gyy[i-1+nx*(j-1+ny*(k-1))] = ( SQ(B_ijk*x_i/exp_nu_ijk)
                                             +SQ(exp_alpha_ijk*y_j) ) /
                                             (SQ(x_i)+SQ(y_j));

                gyz[i-1+nx*(j-1+ny*(k-1))] = 0.0;

                gzz[i-1+nx*(j-1+ny*(k-1))] = SQ(exp_alpha_ijk);


                kxx[i-1+nx*(j-1+ny*(k-1))] = (  ( SQ(r_bar_ijk)*y_j*omega_dx+
                                             (x_i*nu_dy-y_j*nu_dx)*SQ(y_j)
                                             *omega_ijk)*SQ(B_ijk)
                                             +(y_j*B_dx-x_i*B_dy)*omega_ijk
                                             *SQ(y_j)*B_ijk
                                             +(y_j*alpha_dx-x_i*alpha_dy)
                                             *omega_ijk*SQ(x_i*exp_alpha_ijk
                                             *exp_nu_ijk))/(SQ(r_bar_ijk
                                             *exp_nu_ijk)*exp_nu_ijk);
                
                kxy[i-1+nx*(j-1+ny*(k-1))] = ( ( 0.5*SQ(r_bar_ijk)*
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
                
                kxz[i-1+nx*(j-1+ny*(k-1))] = 0.5*SQ(B_ijk)*y_j*omega_dz/
                                             ( SQ(exp_nu_ijk)*exp_nu_ijk );
   
                kyy[i-1+nx*(j-1+ny*(k-1))] = ( ( -SQ(r_bar_ijk)*x_i*omega_dy+
                                             (x_i*nu_dy-y_j*nu_dx)*SQ(x_i)* 
                                             omega_ijk )*SQ(B_ijk) 
                                             +(y_j*B_dx-x_i*B_dy)*omega_ijk
                                             *SQ(x_i)*B_ijk
                                             +(y_j*alpha_dx-x_i*alpha_dy)
                                             *omega_ijk*SQ(y_j*exp_alpha_ijk
                                             *exp_nu_ijk))/(SQ(r_bar_ijk
                                             *exp_nu_ijk)*exp_nu_ijk);
               
                kyz[i-1+nx*(j-1+ny*(k-1))] = -0.5*SQ(B_ijk)*x_i*omega_dz/
                                             ( SQ(exp_nu_ijk)*exp_nu_ijk );
 
                kzz[i-1+nx*(j-1+ny*(k-1))] = (y_j*alpha_dx-x_i*alpha_dy)*
                                             omega_ijk*SQ(exp_alpha_ijk)/
                                             exp_nu_ijk;
              }
 
              press[i-1+nx*(j-1+ny*(k-1))] = pressure_ijk;

              betax[i-1+nx*(j-1+ny*(k-1))] = omega_ijk*y_j;
              betay[i-1+nx*(j-1+ny*(k-1))] = -omega_ijk*x_i;

              betaz[i-1+nx*(j-1+ny*(k-1))] = 0.0;

              eps[i-1+nx*(j-1+ny*(k-1))] = energy_ijk/rho_0_ijk-1.0;
              Y_e[i-1+nx*(j-1+ny*(k-1))] = Y_e_ijk;
              h_ijk = (energy_ijk+pressure_ijk)/rho_0_ijk;

              gamma_ijk = SQ(exp_alpha_ijk)*SQ(B_ijk*exp_alpha_ijk/
                          exp_nu_ijk);

              W_ijk = 1.0/sqrt(1.0-SQ((omega_ijk-Omega_ijk)*B_ijk*
                      r_bar_ijk/SQ(exp_nu_ijk)));


              w_lorentz[i-1+nx*(j-1+ny*(k-1))] = W_ijk;


              velx[i-1+nx*(j-1+ny*(k-1))] = (omega_ijk-Omega_ijk)
                                              *y_j/exp_nu_ijk;

              vely[i-1+nx*(j-1+ny*(k-1))] = -(omega_ijk-Omega_ijk)
                                              *x_i/exp_nu_ijk;

              velz[i-1+nx*(j-1+ny*(k-1))] = 0.0;


              rho[i-1+nx*(j-1+ny*(k-1))] = rho_0_ijk;

                /* if (fabs(x_i - 0.5) < 0.3 && fabs(y_j - 0.5) < 0.3 && fabs(z_k - 0.5) < 0.3) */
                /*   printf("x %g y %g z %g B %g vx %g vy %g alp %g omega %g Omega %g rho %g\n", x_i, y_j, z_k, B_ijk, velx[i-1+nx*(j-1+ny*(k-1))], vely[i-1+nx*(j-1+ny*(k-1))], alp[i-1+nx*(j-1+ny*(k-1))], omega_ijk, Omega_ijk, rho_0_ijk); */



              dens_atm = sqrt(gamma_ijk)*rho_0_atm;
              tau_atm = sqrt(gamma_ijk)*eos_k*pow(rho_0_atm,eos_ideal_fluid_gamma) /
                (eos_ideal_fluid_gamma - 1.0);

              /* *************************************** */
              /* ATMOSPERE SETTINGs                      */
              /* *************************************** */
              if ( (rho[i-1+nx*(j-1+ny*(k-1))] < (1.0 + RNS_atmo_tolerance) * rho_0_atm) ) {
                velx[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                vely[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                velz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                w_lorentz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
              }
              /* *************************************** */
              /* ATMOSPERE SETTINGs                      */
              /* *************************************** */


           } /* END FOR LOOP OF LINE 525 */


      CCTK_INFO("done with loop");
      if( strcmp(zero_shift, "yes")==0) {

        /* SET SHIFT TO ZERO */
 
#pragma omp parallel for
        for(i=1;i<=nx;i++)
           for(j=1;j<=ny;j++)
              for(k=1;k<=nz;k++) {
                 betax[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                 betay[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                 betaz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
              }
      }

      /* compute central value of 3-determinant of metric */

      CCTK_INFO("central value of 3-det");
      *gamma_center = SQ(B[1][1])*exp( 4.0*alpha[1][1] - 2.0*nu[1][1]);

      /* FREE MEMORY */
      
      CCTK_INFO("free mem");
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

      free(s_gp);
      free(mu);



      CCTK_INFO("done");
      return;

} /* END RNSID */


