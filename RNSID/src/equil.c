/****************************************************************************
equil.c		PROCEDURES USED TO COMPUTE THE EQUILIBRIUM STAR
		make_grid(): Sets up the computational grid.
		load_eos(): Loads EOS file.
		e_at_p() etc.: Returns energy density given pressure, etc.
		make_center(): Computes the central pressure and enthalpy
		comp_values(): Given the metric of an equilibrium star,
			computes the mass (M), baryon mass (M_0), kinetic (T)
			 and potential (W) energies, limiting angular velocity
			(Omega_K) and equatorial radius (R_e).
		dm_dr_is() etc.: Isotropic radial derivative of the
			mass function, etc.
		integrate(): Integrate the TOV equations
		guess(): Converts the solution of the TOV equations to the
			right coordinates and uses the spherical solution
			as the first guess for the rotating solution.
		iterate(): The "spin" cycle; iterates the initial guess
			until the field equations are satisfied.
			*****************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "consts.h"
#include "equil.h"
#include "equil_util.h"
#include "rnsid_util.h"

#ifndef RNS_SEQ_COMPILATION
#include <cctk.h>
#else
#define CCTK_VERROR(...) do {fprintf(stderr, __VA_ARGS__); exit(1);} while(0)
#endif

/*******************************************************************/
/* Create computational grid.                                      */
/* Points in the mu-direction are stored in the array mu[i].       */
/* Points in the s-direction are stored in the array s_pg[j].      */
/*******************************************************************/
void make_grid(double s_gp[SDIV+1],
               double mu[MDIV+1])
{
  int m, s;                         /* counters */

      for(s=1;s<=SDIV;s++)
         s_gp[s] = SMAX*(s-1.0)/(SDIV-1.0);

      for(m=1;m<=MDIV;m++)
         mu[m] = (m-1.0)/(MDIV-1.0);
}



/*************************************************************************/
/* Load EOS file.                                                        */
/*************************************************************************/
void load_eos( char eos_file[],
               double log_e_table[MAX_NTAB],
               double log_p_table[MAX_NTAB],
               double log_h_table[MAX_NTAB],
               double log_n0_table[MAX_NTAB],
               int *n_tab)
{
 int i;                    /* counter */

 double p,                 /* pressure */
        rho,               /* density */
        h,                 /* enthalpy */
        n0;                /* number density */

 FILE *f_eos;              /* pointer to eos_file */


    /* OPEN FILE TO READ */

    if((f_eos=fopen(eos_file,"r")) == NULL ) {
       CCTK_VERROR("cannot open file:  %s\n",eos_file);
    }

    /* READ NUMBER OF TABULATED POINTS */

    fscanf(f_eos,"%d\n",n_tab);


    /* READ ENERGY DENSITY, P, H, N0 AND CONVERT TO CACTUS UNITS */

    for(i=1;i<=(*n_tab);i++) {
       fscanf(f_eos,"%lf %lf %lf %lf \n",&rho,&p,&h,&n0) ;
       log_e_table[i]=log10(rho*cactusM*cactusV);     /* multiply by C^2 to get energy density */
       log_p_table[i]=log10(p*cactusM/(cactusL*cactusT*cactusT));
       log_h_table[i]=log10(h*cactusL*cactusL/(cactusT*cactusT));
       log_n0_table[i]=log10(n0*cactusV);
//       printf("WUTANG EOS [e]:%e    [p]:%e    [h]:%e    [n0]:%e \n",log_e_table[i], log_p_table[i], log_h_table[i], log_n0_table[i]);
    }
}
/*************************************************************************/
/* Load Beta equil file.                                                        */
/*************************************************************************/
void load_beta_equil( char beta_equil_file[],
               double log_e_table[MAX_NTAB],
               double Y_e_table[MAX_NTAB],
               int *n_tab_beta)
{
 int i;                    /* counter */

 double rho,               /* density */
        ye;                /* electron fraction */

 FILE *f_beta;              /* pointer to beta_equil_file */


    /* OPEN FILE TO READ */

    if((f_beta=fopen(beta_equil_file,"r")) == NULL ) {
       CCTK_VERROR("cannot open beta-equil. file:  %s\n",beta_equil_file);
    }

    /* READ NUMBER OF TABULATED POINTS */

    fscanf(f_beta,"%d\n",n_tab_beta);

    /* READ ENERGY DENSITY, P, H, N0 AND CONVERT TO CACTUS UNITS */

    for(i=1;i<=(*n_tab_beta);i++) {
       fscanf(f_beta,"%lf %lf \n",&rho,&ye) ;
       log_e_table[i]=log10(rho*cactusM*cactusV);     /* multiply by C^2 to get energy density */
       Y_e_table[i]=(ye);
//       printf("WUTANG EOS [e]:%e    [p]:%e    [h]:%e    [n0]:%e \n",log_e_table[i], log_p_table[i], log_h_table[i], log_n0_table[i]);
    }
}



/*******************************************************************/
double e_of_rho0(double rho0, double Gamma_P, double eos_k)
{
 return(eos_k*pow(rho0,Gamma_P)/(Gamma_P-1.0)+rho0);
}


/*C*/
/*******************************************************************/
double e_at_p(double pp,
              double log_e_tab[MAX_NTAB],
              double log_p_tab[MAX_NTAB],
              int    n_tab,
              int    *n_nearest_pt,
              char eos_type[],
	      double eos_k,
              double Gamma_P)
{
 if(strcmp(eos_type,"tab")==0)
   return pow(10.0,interp(log_p_tab,log_e_tab,n_tab,log10(pp), n_nearest_pt));
 else
   return pp/(Gamma_P-1.0) + pow(pp/eos_k,1.0/Gamma_P);
}

/*C*/
/*******************************************************************/
double p_at_e(double ee,
              double log_p_tab[MAX_NTAB],
              double log_e_tab[MAX_NTAB],
              int    n_tab,
              int    *n_nearest_pt)
{
 return pow(10.0,interp(log_e_tab,log_p_tab,n_tab,log10(ee), n_nearest_pt));
}

/*C*/
/*******************************************************************/
double p_at_h(double hh,
              double log_p_tab[MAX_NTAB],
              double log_h_tab[MAX_NTAB],
              int    n_tab,
              int    *n_nearest_pt)
{
 return pow(10.0,interp(log_h_tab,log_p_tab,n_tab,log10(hh), n_nearest_pt));
}

/*C*/
/*******************************************************************/
double h_at_p(double pp,
              double log_h_tab[MAX_NTAB],
              double log_p_tab[MAX_NTAB],
              int    n_tab,
              int    *n_nearest_pt)
{
 return pow(10.0,interp(log_p_tab,log_h_tab,n_tab,log10(pp), n_nearest_pt));
}

/*C*/
/*******************************************************************/
double n0_at_e(double ee,
               double log_n0_tab[MAX_NTAB],
               double log_e_tab[MAX_NTAB],
               int    n_tab,
               int    *n_nearest_pt)
{
 return pow(10.0,interp(log_e_tab,log_n0_tab,n_tab,log10(ee), n_nearest_pt));
}
/*C*/
/***************************************************************/
double e_at_n0(double nn0,
               double log_e_tab[MAX_NTAB],
               double log_n0_tab[MAX_NTAB],
               int    n_tab,
               int    *n_nearest_pt)
{
 return pow(10.0,interp(log_n0_tab,log_e_tab,n_tab,log10(nn0), n_nearest_pt));
}
/*C*/
/***************************************************************/
void make_center(double e_center,
                 double log_e_tab[MAX_NTAB],
                 double log_p_tab[MAX_NTAB],
                 double log_h_tab[MAX_NTAB],
                 int n_tab,
                 char eos_type[],
		 double eos_k,
                 double Gamma_P,
                 double *p_center,
                 double *h_center)
{
 int n_nearest;

 double rho0_center;

 n_nearest=n_tab/2;

 if(strcmp(eos_type,"tab")==0) {
   (*p_center) = p_at_e( e_center, log_p_tab, log_e_tab, n_tab, &n_nearest);
   (*h_center) = h_at_p( (*p_center), log_h_tab, log_p_tab, n_tab, &n_nearest);
 }
 else {
       rho0_center = rtsec_G( e_of_rho0, Gamma_P, 0.0,e_center,DBL_EPSILON,
                              e_center, eos_k );
       (*p_center) = eos_k*pow(rho0_center,Gamma_P);
       (*h_center) = log((e_center+(*p_center))/rho0_center);
 }


}

/*C*/
/********************************************************************/
/* Compute various quantities.                                      */
/********************************************************************/
void comp_values(double s_gp[SDIV+1],
                 double mu[MDIV+1],
                 double r_ratio,
                 double e_surface,
                 double r_e,
                 char   eos_type[],
                 double log_e_tab[MAX_NTAB],
                 double log_n0_tab[MAX_NTAB],
                 int    n_tab,
                 double Omega,
                 double **rho,
                 double **gama,
                 double **alpha,
                 double **omega,
                 double **energy,
                 double **pressure,
                 double **enthalpy,
                 double **velocity_sq,
                 double *Mass,
                 double *Mass_0,
                 double *T,
                 double *W,
                 double *Omega_K,
                 double *R_e,
                 char   rotation_type[],
		 double **Omega_diff,
		 double *J)
{
 int s,
     m,
     n_nearest;


 double **velocity,                /* velocity */
        **rho_0,                   /* rest mass density */
        gama_equator,              /* gama at equator */
        rho_equator,               /* rho at equator */
        omega_equator,             /* omega at equator */
        D_m[SDIV+1],               /* int. quantity for M */
        D_m_0[SDIV+1],             /* int. quantity for M_0 */
        D_m_p[SDIV+1],             /* int. quantity for M_p */
        D_J[SDIV+1],               /* int. quantity for J */
        d_o_e[SDIV+1],
        d_g_e[SDIV+1],
        d_r_e[SDIV+1],
        doe,
        dge,
        dre,
        vek,
        s_e,
        gama_mu_0[SDIV+1],
        rho_mu_0[SDIV+1],
        omega_mu_0[SDIV+1],
        Mass_p,
        T_diff,
        D_T[SDIV+1];               /* int. quantity for J */


   s_e=0.5;

   velocity = array_allocate(1,SDIV,1,MDIV);
   rho_0 = array_allocate(1,SDIV,1,MDIV);

   for(s=1;s<=SDIV;s++)
    for(m=1;m<=MDIV;m++)
     velocity[s][m]=sqrt(velocity_sq[s][m]);

   for(s=1;s<=SDIV;s++) {
      gama_mu_0[s]=gama[s][1];
      rho_mu_0[s]=rho[s][1];
      omega_mu_0[s]=omega[s][1];
   }

   n_nearest= SDIV/2;
   gama_equator=interp(s_gp,gama_mu_0,SDIV,s_e, &n_nearest);
   rho_equator=interp(s_gp,rho_mu_0,SDIV,s_e, &n_nearest);

   if(r_ratio==1.0) {
    omega_equator=0.0;
   }
   else {
    omega_equator=interp(s_gp,omega_mu_0,SDIV,s_e, &n_nearest);
   }

/* Circumferential radius */

   if(strcmp(eos_type,"tab")==0)
       (*R_e) = r_e*exp((gama_equator-rho_equator)/2.0);
/*       (*R_e) = sqrt(KAPPA)*r_e*exp((gama_equator-rho_equator)/2.0); */ /*conversion*/
   else
     (*R_e) = r_e*exp((gama_equator-rho_equator)/2.0);

 /* Masses and angular momentum */

   (*Mass) = 0.0;              /* initialize */
   (*Mass_0) = 0.0;
   Mass_p = 0.0;
   (*J)=0.0;
   T_diff = 0.0;

 if(strcmp(eos_type,"tab")==0) {
   n_nearest=n_tab/2;
   for(s=1;s<=SDIV;s++)
      for(m=1;m<=MDIV;m++) {
           if(energy[s][m]>e_surface)
/*             rho_0[s][m]=n0_at_e(energy[s][m], log_n0_tab, log_e_tab, n_tab,
                                             &n_nearest)*MB*KSCALE*SQ(C); */
             /*this conversion factor should be MB*cactusM*/
	     rho_0[s][m]=n0_at_e(energy[s][m], log_n0_tab, log_e_tab, n_tab,
				 &n_nearest)*MB*cactusM;
           else
             rho_0[s][m]=0.0;
      }
 }else {
        for(s=1;s<=SDIV;s++)
           for(m=1;m<=MDIV;m++)
                rho_0[s][m]=(energy[s][m]+pressure[s][m])*exp(-enthalpy[s][m]);
  }


   for(s=1;s<=SDIV;s++) {
    D_m[s]=0.0;           /* initialize */
    D_m_0[s]=0.0;
    D_m_p[s]=0.0;
    D_J[s]=0.0;
    D_T[s]=0.0;

    for(m=1;m<=MDIV-2;m+=2) {
     D_m[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+gama[s][m])*
              (((energy[s][m]+pressure[s][m])/(1.0-velocity_sq[s][m]))*
              (1.0+velocity_sq[s][m]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m]*mu[m])*r_e*omega[s][m]*
              exp(-rho[s][m])) + 2.0*pressure[s][m])

            + 4.0*exp(2.0*alpha[s][m+1]+gama[s][m+1])*
              (((energy[s][m+1]+pressure[s][m+1])/(1.0-velocity_sq[s][m+1]))*
              (1.0+velocity_sq[s][m+1]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+1])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+1]*mu[m+1])*r_e*omega[s][m+1]*
              exp(-rho[s][m+1])) + 2.0*pressure[s][m+1])

            + exp(2.0*alpha[s][m+2]+gama[s][m+2])*
              (((energy[s][m+2]+pressure[s][m+2])/(1.0-velocity_sq[s][m+2]))*
              (1.0+velocity_sq[s][m+2]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+2])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+2]*mu[m+2])*r_e*omega[s][m+2]*
              exp(-rho[s][m+2])) + 2.0*pressure[s][m+2]));

     D_m_0[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+(gama[s][m]
              -rho[s][m])/2.0)*rho_0[s][m]/sqrt(1.0-velocity_sq[s][m])

             + 4.0* exp(2.0*alpha[s][m+1]+(gama[s][m+1]
             -rho[s][m+1])/2.0)*rho_0[s][m+1]/sqrt(1.0-velocity_sq[s][m+1])

             + exp(2.0*alpha[s][m+2]+(gama[s][m+2]
             -rho[s][m+2])/2.0)*rho_0[s][m+2]/sqrt(1.0-velocity_sq[s][m+2]));

     D_m_p[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+(gama[s][m]
              -rho[s][m])/2.0)*energy[s][m]/sqrt(1.0-velocity_sq[s][m])

              + 4.0* exp(2.0*alpha[s][m+1]+(gama[s][m+1]
              -rho[s][m+1])/2.0)*energy[s][m+1]/sqrt(1.0-velocity_sq[s][m+1])

             + exp(2.0*alpha[s][m+2]+(gama[s][m+2]
             -rho[s][m+2])/2.0)*energy[s][m+2]/sqrt(1.0-velocity_sq[s][m+2]));

     D_J[s] += (1.0/(3.0*(MDIV-1)))*( sqrt(1.0-mu[m]*mu[m])*
              exp(2.0*alpha[s][m]+gama[s][m]-rho[s][m])*(energy[s][m]
              +pressure[s][m])*sqrt(velocity_sq[s][m])/(1.0-velocity_sq[s][m])

              +4.0*sqrt(1.0-mu[m+1]*mu[m+1])*
              exp(2.0*alpha[s][m+1]+gama[s][m+1]-rho[s][m+1])*(energy[s][m+1]
              +pressure[s][m+1])*sqrt(velocity_sq[s][m+1])/
              (1.0-velocity_sq[s][m+1])

              + sqrt(1.0-mu[m+2]*mu[m+2])*
              exp(2.0*alpha[s][m+2]+gama[s][m+2]-rho[s][m+2])*(energy[s][m+2]
              +pressure[s][m+2])*sqrt(velocity_sq[s][m+2])/
              (1.0-velocity_sq[s][m+2]));

     D_T[s] += (1.0/(3.0*(MDIV-1)))*( sqrt(1.0-mu[m]*mu[m])*
              exp(2.0*alpha[s][m]+gama[s][m]-rho[s][m])*(energy[s][m]
              +pressure[s][m])*sqrt(velocity_sq[s][m])*Omega_diff[s][m]/
              (1.0-velocity_sq[s][m])

              +4.0*sqrt(1.0-mu[m+1]*mu[m+1])*
              exp(2.0*alpha[s][m+1]+gama[s][m+1]-rho[s][m+1])*(energy[s][m+1]
              +pressure[s][m+1])*sqrt(velocity_sq[s][m+1])*Omega_diff[s][m+1]/
              (1.0-velocity_sq[s][m+1])

              + sqrt(1.0-mu[m+2]*mu[m+2])*
              exp(2.0*alpha[s][m+2]+gama[s][m+2]-rho[s][m+2])*(energy[s][m+2]
              +pressure[s][m+2])*sqrt(velocity_sq[s][m+2])*Omega_diff[s][m+2]/
              (1.0-velocity_sq[s][m+2]));
    }
   }

    for(s=1;s<=SDIV-4;s+=2) {
     (*Mass) += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_m[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m[s+2]);

     (*Mass_0) += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_m_0[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m_0[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m_0[s+2]);

     Mass_p += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_m_p[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m_p[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m_p[s+2]);

     (*J) += (SMAX/(3.0*(SDIV-1)))*((pow(s_gp[s],3.0)/pow(1.0-s_gp[s],5.0))*
          D_J[s]+ 4.0*(pow(s_gp[s+1],3.0)/pow(1.0-s_gp[s+1],5.0))*
          D_J[s+1] + (pow(s_gp[s+2],3.0)/pow(1.0-s_gp[s+2],5.0))*
          D_J[s+2]);

     T_diff += (SMAX/(3.0*(SDIV-1)))*((pow(s_gp[s],3.0)/pow(1.0-s_gp[s],5.0))*
          D_T[s]+ 4.0*(pow(s_gp[s+1],3.0)/pow(1.0-s_gp[s+1],5.0))*
          D_T[s+1] + (pow(s_gp[s+2],3.0)/pow(1.0-s_gp[s+2],5.0))*
          D_T[s+2]);
    }

    if(strcmp(eos_type,"tab")==0) {
/*	(*Mass) *= 4*PI*sqrt(KAPPA)*C*C*pow(r_e,3.0)/G; conversions
	(*Mass_0) *= 4*PI*sqrt(KAPPA)*C*C*pow(r_e,3.0)/G;
	Mass_p *= 4*PI*sqrt(KAPPA)*C*C*pow(r_e,3.0)/G;*/
	(*Mass) *= 4*PI*pow(r_e,3.0);
	(*Mass_0) *= 4*PI*pow(r_e,3.0);
	Mass_p *= 4*PI*pow(r_e,3.0);
    }
    else {
          (*Mass) *= 4*PI*pow(r_e,3.0);
          (*Mass_0) *= 4*PI*pow(r_e,3.0);
          Mass_p *= 4*PI*pow(r_e,3.0);
    }

    if(r_ratio==1.0)
         (*J)=0.0;
    else {
          if(strcmp(eos_type,"tab")==0) {
/*              (*J) *= 4*PI*KAPPA*C*C*C*pow(r_e,4.0)/G;  conversions
              T_diff *= 2.0*PI*sqrt(KAPPA)*C*C*C*C*pow(r_e,4.0)/G;*/
              (*J) *= 4*PI*pow(r_e,4.0);
              T_diff *= 2.0*PI*pow(r_e,4.0);
          }
          else {
              (*J) *= 4*PI*pow(r_e,4.0);
              T_diff *= 2.0*PI*pow(r_e,4.0);
          }
    }
    if(strcmp(rotation_type,"uniform")==0)
        (*T)=0.5*(*J)*Omega;
    else
        (*T)= T_diff;


    if(strcmp(eos_type,"tab")==0)
/*	(*W) = Mass_p*C*C - (*Mass)*C*C + (*T); */ /* conversions*/
	(*W) = Mass_p - (*Mass) + (*T);
    else
      (*W) = Mass_p - (*Mass) + (*T);


/* Kepler angular velocity */

/* MODIFY FOR DIFF. ROTATION */

   for(s=1;s<=SDIV;s++) {
     d_o_e[s]=deriv_s(omega,s,1);
     d_g_e[s]=deriv_s(gama,s,1);
     d_r_e[s]=deriv_s(rho,s,1);
   }

   n_nearest=SDIV/2;
   doe=interp(s_gp,d_o_e,SDIV,s_e, &n_nearest);
   dge=interp(s_gp,d_g_e,SDIV,s_e, &n_nearest);
   dre=interp(s_gp,d_r_e,SDIV,s_e, &n_nearest);

  vek=(doe/(8.0+dge-dre))*r_e*exp(-rho_equator) + sqrt(((dge+dre)/(8.0+dge
        -dre)) + pow((doe/(8.0+dge-dre))*r_e*exp(-rho_equator),2.0));

  if(strcmp(eos_type,"tab")==0)
/*      (*Omega_K) = (C/sqrt(KAPPA))*(omega_equator+vek*exp(rho_equator)/r_e); */ /*conversions*/
      (*Omega_K) = (omega_equator+vek*exp(rho_equator)/r_e);
  else
    (*Omega_K) = omega_equator + vek*exp(rho_equator)/r_e;


   array_free(velocity,1,SDIV,1,MDIV);
   array_free(rho_0,1,SDIV,1,MDIV);

}

/*C*/
/**************************************************************************/
double dm_dr_is(double r_is,
                double r,
                double m,
                double p,
                double e_center,
                double p_surface,
                double log_e_tab[MAX_NTAB],   /* was [SDIV+1], corrected by niksterg on 16.10.14 */
                double log_p_tab[MAX_NTAB],   /* was [SDIV+1], corrected by niksterg on 16.10.14 */
                int    n_tab,
                int    *n_nearest_pt,
                char eos_type[],
		double eos_k,
                double Gamma_P)
{
 double dmdr,
        e_d;

 if(p<p_surface)
    e_d=0.0;
 else
   e_d = e_at_p(p, log_e_tab, log_p_tab, n_tab, n_nearest_pt, eos_type, eos_k, Gamma_P);

 if(r_is<RMIN)
    dmdr=4.0*PI*e_center*r*r*(1.0+4.0*PI*e_center*r*r/3.0);
 else
    dmdr=4.0*PI*e_d*r*r*r*sqrt(1.0-2.0*m/r)/r_is;

return dmdr;
}

/*C*/
/**************************************************************************/
double dp_dr_is(double r_is,
                double r,
                double m,
                double p,
                double e_center,
                double p_surface,
                double log_e_tab[MAX_NTAB],   /* was [SDIV+1], corrected by niksterg on 160016.10.14 */
                double log_p_tab[MAX_NTAB],   /* was [SDIV+1], corrected by niksterg on 160016.10.14 */
                int    n_tab,
                int    *n_nearest_pt,
                char eos_type[],
		double eos_k,
                double Gamma_P)
{ double dpdr,
         e_d;

  if(p<p_surface) e_d=0.0;
  else
    e_d=e_at_p(p, log_e_tab, log_p_tab, n_tab, n_nearest_pt, eos_type, eos_k, Gamma_P);

  if(r_is<RMIN) dpdr = -4.0*PI*(e_center+p)*(e_center+3.0*p)*r*(1.0
                     +4.0*e_center*r*r/3.0)/3.0;

  else
   dpdr = -(e_d+p)*(m+4.0*PI*r*r*r*p)/(r*r_is*sqrt(1.0-2.0*m/r));

 return dpdr;
}


/**************************************************************************/
double dr_dr_is(double r_is, double r, double m)
{
 double drdris;

 if(r_is<RMIN) drdris=1.0;
  else
   drdris=(r/r_is)*sqrt(1.0-2.0*m/r);

 return drdris;
}

/*C*/
/************************************************************************/
void integrate(int    i_check,
               char   eos_type[],
	       double eos_k,
               double e_center,
               double p_center,
               double p_surface,
               double e_surface,
               double Gamma_P,
               double log_e_tab[MAX_NTAB],
               double log_p_tab[MAX_NTAB],
               double log_h_tab[MAX_NTAB],
               int    n_tab,
               double r_is_gp[RDIV+1],
               double lambda_gp[RDIV+1],
               double nu_gp[RDIV+1],
               double *r_is_final,
               double *r_final,
               double *m_final)
{
  int i=2,
      n_nearest;

  double r,                           /* radius */
         r_is,                        /* isotropic radial coordinate */
         r_is_est,                    /* estimate on final isotr. radius */
         r_is_check,                  /*                      */
         dr_is_save,                  /* r_is saving interval */
         rho0,
         e_d,                         /* density */
         p,                           /* pressure */
         h,                           /* stepsize during integration */
         m,                           /* mass   */
         nu_s,
         hh,
         a1,a2,a3,a4,b1,b2,b3,b4,     /* coeff. in Runge-Kutta equations */
         c1,c2,c3,c4,
         k_rescale,
         r_gp[RDIV+1],
         e_d_gp[RDIV+1];


    if(i_check==1) {
      if(strcmp(eos_type,"tab")==0)
	  r_is_est=1.5e6/sqrt(KAPPA);/*is this a conversion or not?*/
      else
        r_is_est=2.0*sqrt(Gamma_P/(4.0*PI*(Gamma_P-1.0)))*
                            pow(e_center,(Gamma_P-2.0)/2.0);

      h=r_is_est/100;
    }
    else {
          r_is_est= (*r_is_final);
          h=r_is_est/10000;
	 }

    dr_is_save = (*r_is_final)/RDIV;
    r_is_check = dr_is_save;
    r_is=0.0;                            /* initial isotropic radius */
    r=0.0;                               /* initial radius */
    m=0.0;                               /* initial mass */
    p=p_center;                          /* initial pressure */

    r_is_gp[1]=0.0;
    r_gp[1]=0.0;
    lambda_gp[1]=0.0;
    e_d_gp[1] = e_center;

    n_nearest = n_tab/2;

    while (p>=p_surface) {

      e_d = e_at_p(p, log_e_tab, log_p_tab, n_tab, &n_nearest, eos_type, eos_k, Gamma_P);

     if((i_check==3) && (r_is>r_is_check) && (i<=RDIV)) {
      r_is_gp[i]=r_is;
      r_gp[i]=r;
      e_d_gp[i]=e_d;
      i++;
      r_is_check += dr_is_save;
     }

     (*r_is_final)=r_is;
     (*r_final)=r;
     (*m_final)=m;



     a1=dr_dr_is(r_is,r,m);

     b1=dm_dr_is(r_is,r,m,p, e_center, p_surface, log_e_tab, log_p_tab, n_tab,
                                                &n_nearest, eos_type, eos_k, Gamma_P);
     c1=dp_dr_is(r_is,r,m,p, e_center, p_surface, log_e_tab, log_p_tab, n_tab,
                                                &n_nearest, eos_type, eos_k, Gamma_P);



     a2=dr_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0);

     b2=dm_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0, p+h*c1/2.0, e_center,
                          p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest,
                          eos_type, eos_k, Gamma_P);

     c2=dp_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0, p+h*c1/2.0, e_center,
                          p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest,
                          eos_type, eos_k, Gamma_P);



     a3=dr_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0);

     b3=dm_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0, p+h*c2/2.0, e_center,
                          p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest,
                          eos_type, eos_k, Gamma_P);

     c3=dp_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0, p+h*c2/2.0, e_center,
                          p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest,
                          eos_type, eos_k, Gamma_P);



     a4=dr_dr_is(r_is+h, r+h*a3, m+h*b3);

     b4=dm_dr_is(r_is+h, r+h*a3, m+h*b3, p+h*c3, e_center, p_surface,
                                log_e_tab, log_p_tab, n_tab,&n_nearest,
                                eos_type, eos_k, Gamma_P);

     c4=dp_dr_is(r_is+h, r+h*a3, m+h*b3, p+h*c3, e_center, p_surface,
                                log_e_tab, log_p_tab, n_tab,&n_nearest,
                                eos_type, eos_k, Gamma_P);



     r += (h/6.0)*(a1+2*a2+2*a3+a4);
     m += (h/6.0)*(b1+2*b2+2*b3+b4);
     p += (h/6.0)*(c1+2*c2+2*c3+c4);

     r_is += h;

    }

    r_is_gp[RDIV]=(*r_is_final);
    r_gp[RDIV]=(*r_final);


/* Rescale r_is and compute lambda */

   if(i_check==3) {
      k_rescale=0.5*((*r_final)/(*r_is_final))*(1.0-(*m_final)/(*r_final)+
                sqrt(1.0-2.0*(*m_final)/(*r_final)));

      (*r_is_final) *= k_rescale;

      nu_s = log((1.0-(*m_final)/(2.0*(*r_is_final)))/(1.0+(*m_final)/
               (2.0*(*r_is_final))));

      for(i=1;i<=RDIV;i++) {
         r_is_gp[i] *= k_rescale;

         if(i==1) lambda_gp[1]= log(1.0/k_rescale);
           else
               lambda_gp[i]=log(r_gp[i]/r_is_gp[i]);

         if(e_d_gp[i]<e_surface) 
 
           hh=0.0;

         else { if(strcmp(eos_type,"tab")==0) {
                  p=p_at_e(e_d_gp[i], log_p_tab, log_e_tab, n_tab, &n_nearest);
                  hh=h_at_p(p, log_h_tab, log_p_tab, n_tab, &n_nearest);
                }
                else {
                     rho0=rtsec_G(e_of_rho0,Gamma_P,0.0,e_d_gp[i],DBL_EPSILON,
				  e_d_gp[i], eos_k);
                      p=eos_k*pow(rho0,Gamma_P);
                      hh=log((e_d_gp[i]+p)/rho0);
                }
              }

         nu_gp[i]=nu_s-hh;
      }
      nu_gp[RDIV]=nu_s;
   }
}

/*C*/
/*************************************************************************/
void guess(double s_gp[SDIV+1],
           char   eos_type[],
	   double eos_k,
           double e_center,
           double p_center,
           double p_surface,
           double e_surface,
           double Gamma_P,
           double log_e_tab[MAX_NTAB],
           double log_p_tab[MAX_NTAB],
           double log_h_tab[MAX_NTAB],
           int    n_tab,
           double **rho,
           double **gama,
           double **alpha,
           double **omega,
           double *r_e)
{
 int s,
     m,
     n_nearest;

 double r_is_s,
        r_is_final=0.0,
        r_final,
        m_final,
        lambda_s,
        nu_s,
        r_is_gp[RDIV+1],
        lambda_gp[RDIV+1],
        nu_gp[RDIV+1],
        gama_mu_0[SDIV+1],
        rho_mu_0[SDIV+1],
        gama_eq,
        rho_eq,
        s_e=0.5;


 integrate(1, eos_type, eos_k, e_center, p_center, p_surface, e_surface, Gamma_P,
              log_e_tab, log_p_tab, log_h_tab, n_tab, r_is_gp, lambda_gp,
              nu_gp, &r_is_final, &r_final, &m_final);

 integrate(2, eos_type, eos_k, e_center, p_center, p_surface, e_surface, Gamma_P,
              log_e_tab, log_p_tab, log_h_tab, n_tab, r_is_gp, lambda_gp,
              nu_gp, &r_is_final, &r_final, &m_final);

 integrate(3, eos_type, eos_k, e_center, p_center, p_surface, e_surface, Gamma_P,
              log_e_tab, log_p_tab, log_h_tab, n_tab, r_is_gp, lambda_gp,
              nu_gp, &r_is_final, &r_final, &m_final);


 n_nearest=RDIV/2;
 for(s=1;s<=SDIV;s++) {
    r_is_s=r_is_final*(s_gp[s]/(1.0-s_gp[s]));

    if(r_is_s<r_is_final) {
      lambda_s=interp(r_is_gp,lambda_gp,RDIV,r_is_s,&n_nearest);
      nu_s=interp(r_is_gp,nu_gp,RDIV,r_is_s,&n_nearest);
    }
    else {
      lambda_s=2.0*log(1.0+m_final/(2.0*r_is_s));
      nu_s=log((1.0-m_final/(2.0*r_is_s))/(1.0+m_final/(2*r_is_s)));
    }

    gama[s][1]=nu_s+lambda_s;
    rho[s][1]=nu_s-lambda_s;

    for(m=1;m<=MDIV;m++) {
        gama[s][m]=gama[s][1];
        rho[s][m]=rho[s][1];
        alpha[s][m]=(gama[s][1]-rho[s][1])/2.0;
        omega[s][m]=0.0;
    }

    gama_mu_0[s]=gama[s][1];                   /* gama at \mu=0 */
    rho_mu_0[s]=rho[s][1];                     /* rho at \mu=0 */

 }

   n_nearest=SDIV/2;
   gama_eq = interp(s_gp,gama_mu_0,SDIV,s_e,&n_nearest); /* gama at equator */
   rho_eq = interp(s_gp,rho_mu_0,SDIV,s_e,&n_nearest);   /* rho at equator */

   (*r_e)= r_final*exp(0.5*(rho_eq-gama_eq));

}


/*C*/
/*************************************************************************/
/* Main iteration cycle.                                                 */
/*************************************************************************/
void iterate(double s_gp[SDIV+1],
             double mu[MDIV+1],
             char   eos_type[],
	     double eos_k,
             double log_e_tab[MAX_NTAB],
             double log_p_tab[MAX_NTAB],
             double log_h_tab[MAX_NTAB],
             int    n_tab,
             double Gamma_P,
             double r_ratio,
             double h_center,
             double enthalpy_min,
             int    a_check,
             double accuracy,
             int    print_dif,
             double cf,
             double *r_e_new,
             double **rho,
             double **gama,
             double **alpha,
             double **omega,
             double **energy,
             double **pressure,
             double **enthalpy,
             double **velocity_sq,
             double *Omega,
             char   rotation_type[],
             double A_diff,
             double *Omega_e,
             double **Omega_diff,
             int RNS_lmax)
 {
 int m,                      /* counter */
     s,                      /* counter */
     n,                      /* counter */
     k,                      /* counter */
     n_of_it=0,              /* number of iterations */
     n_nearest,
     i,
     j;

 /* ------------------------------
    iteration counter: it stop
    when  MAX_ITERATION  = 200;
    is reached. It is definined
    in const.h
    ---------------------------- */
int     ITERATION_COUNT= 0;

double  ***f_rho,
        ***f_gama;


double   sum_rho=0.0,         /* intermediate sum in eqn for rho */
	 sum_gama=0.0,        /* intermediate sum in eqn for gama */
	 sum_omega=0.0,       /* intermediate sum in eqn for omega */
         r_e_old,             /* equatorial radius in previus cycle */
   	 dif,                 /* difference | r_e_old - r_e | */
         d_gama_s,            /* derivative of gama w.r.t. s */
         d_gama_m,            /* derivative of gama w.r.t. m */
         d_rho_s,             /* derivative of rho w.r.t. s */
         d_rho_m,             /* derivative of rho w.r.t. m */
         d_omega_s,           /* derivative of omega w.r.t. s */
         d_omega_m,           /* derivative of omega w.r.t. m */
         d_gama_ss,           /* 2nd derivative of gama w.r.t. s */
         d_gama_mm,           /* 2nd derivative of gama w.r.t. m */
         d_gama_sm,           /* derivative of gama w.r.t. m and s */
         temp1,                /* temporary term in da_dm */
         temp2,
         temp3,
         temp4,
         temp5,
         temp6,
         temp7,
         temp8,
         m1,
         s1,
         s2,
         ea,
         rsm,
         gsm,
         omsm,
         esm,
         psm,
         v2sm,
         mum,
         sgp,
         s_1,
         e_gsm,
         e_rsm,
         rho0sm,
         term_in_Omega_h,
         r_p,
         s_p,
         gama_pole_h,                  /* gama^hat at pole */
         gama_center_h,                /* gama^hat at center */
         gama_equator_h,               /* gama^hat at equator */
         rho_pole_h,                   /* rho^hat at pole */
         rho_center_h,                 /* rho^hat at center */
         rho_equator_h,                /* rho^hat at equator */
         omega_equator_h,              /* omega^hat at equator */
         gama_mu_1[SDIV+1],            /* gama at \mu=1 */
         gama_mu_0[SDIV+1],            /* gama at \mu=0 */
         rho_mu_1[SDIV+1],             /* rho at \mu=1 */
         rho_mu_0[SDIV+1],             /* rho at \mu=0 */
         omega_mu_0[SDIV+1],           /* omega at \mu=0 */
         s_e=0.5,
       **da_dm,
       **dgds,
       **dgdm,
       **D1_rho,
       **D1_gama,
       **D1_omega,
       **D2_rho,
       **D2_gama,
       **D2_omega,
       **S_gama,
       **S_rho,
       **S_omega,
       **f2n,
       **P_2n,
       **P1_2n_1,
         Omega_h,
         sin_theta[MDIV+1],
         theta[MDIV+1],
         sk,
         sj,
         sk1,
         sj1,
         r_e,
         d_e_ds_mu_0[SDIV/2],
         s_max,
         Omega_lower_bracket,
         Omega_upper_bracket,
         Omega_e_lower_bracket,
         Omega_e_upper_bracket,
         Omega_c,
         tol;

    const int is_uniform_rotation = strcmp(rotation_type,"uniform")==0;

    f2n = array_allocate(1,RNS_lmax+1,1,SDIV);
    f_rho = tensor_allocate(1,SDIV,1,RNS_lmax+1,1,SDIV);
    f_gama = tensor_allocate(1,SDIV,1,RNS_lmax+1,1,SDIV);

    P_2n = array_allocate(1,MDIV,1,RNS_lmax+1);
    P1_2n_1 = array_allocate(1,MDIV,1,RNS_lmax+1);


    for(n=0;n<=RNS_lmax;n++)
       for(i=2;i<=SDIV;i++) f2n[n+1][i] = pow((1.0-s_gp[i])/s_gp[i],2.0*n);

    if(SMAX!=1.0) {

     for(j=2;j<=SDIV;j++)
        for(n=1;n<=RNS_lmax;n++)
           for(k=2;k<=SDIV;k++) {
                 sk=s_gp[k];
                 sj=s_gp[j];
                 sk1=1.0-sk;
                 sj1=1.0-sj;

                 if(k<j) {
                          f_rho[j][n+1][k] = f2n[n+1][j]*sj1/(sj*
                                  f2n[n+1][k]*sk1*sk1);
                          f_gama[j][n+1][k] = f2n[n+1][j]/(f2n[n+1][k]*sk*sk1);
	                 }else {
                          f_rho[j][n+1][k] = f2n[n+1][k]/(f2n[n+1][j]*sk*sk1);
                          f_gama[j][n+1][k] = f2n[n+1][k]*sj1*sj1*sk/(sj*sj
                                            *f2n[n+1][j]*sk1*sk1*sk1);
                 }
	    }
     j=1;

       n=0;
       for(k=2;k<=SDIV;k++) {
          sk=s_gp[k];
          f_rho[j][n+1][k]=1.0/(sk*(1.0-sk));
       }

       n=1;
       for(k=2;k<=SDIV;k++) {
          sk=s_gp[k];
          sk1=1.0-sk;
          f_rho[j][n+1][k]=0.0;
          f_gama[j][n+1][k]=1.0/(sk*sk1);
       }

       for(n=2;n<=RNS_lmax;n++)
          for(k=1;k<=SDIV;k++) {
             f_rho[j][n+1][k]=0.0;
             f_gama[j][n+1][k]=0.0;
          }


     k=1;

       n=0;
       for(j=1;j<=SDIV;j++)
          f_rho[j][n+1][k]=0.0;

       for(j=1;j<=SDIV;j++)
          for(n=1;n<=RNS_lmax;n++) {
             f_rho[j][n+1][k]=0.0;
             f_gama[j][n+1][k]=0.0;
          }


     n=0;
     for(j=2;j<=SDIV;j++)
        for(k=2;k<=SDIV;k++) {
               sk=s_gp[k];
               sj=s_gp[j];
               sk1=1.0-sk;
               sj1=1.0-sj;

               if(k<j)
                 f_rho[j][n+1][k] = sj1/(sj*sk1*sk1);
               else
                 f_rho[j][n+1][k] = 1.0/(sk*sk1);
      }

   }
   else{      
        for(j=2;j<=SDIV-1;j++)
           for(n=1;n<=RNS_lmax;n++)
              for(k=2;k<=SDIV-1;k++) {
                 sk=s_gp[k];
                 sj=s_gp[j];
                 sk1=1.0-sk;
                 sj1=1.0-sj;

                 if(k<j) {   
                          f_rho[j][n+1][k] = f2n[n+1][j]*sj1/(sj*
                                           f2n[n+1][k]*sk1*sk1);
                          f_gama[j][n+1][k] = f2n[n+1][j]/(f2n[n+1][k]*sk*sk1);
                 }else {     
                          f_rho[j][n+1][k] = f2n[n+1][k]/(f2n[n+1][j]*sk*sk1);

                          f_gama[j][n+1][k] = f2n[n+1][k]*sj1*sj1*sk/(sj*sj
                                            *f2n[n+1][j]*sk1*sk1*sk1);
                 }
	      }
   
        j=1;
 
          n=0; 
          for(k=2;k<=SDIV-1;k++) {
             sk=s_gp[k];
             f_rho[j][n+1][k]=1.0/(sk*(1.0-sk));
          }

          n=1;
          for(k=2;k<=SDIV-1;k++) {
             sk=s_gp[k];
             sk1=1.0-sk;         
             f_rho[j][n+1][k]=0.0;
             f_gama[j][n+1][k]=1.0/(sk*sk1);
          }

          for(n=2;n<=RNS_lmax;n++)
             for(k=1;k<=SDIV-1;k++) {
                f_rho[j][n+1][k]=0.0;
                f_gama[j][n+1][k]=0.0;
             }

        k=1;
 
          n=0;
          for(j=1;j<=SDIV-1;j++)
             f_rho[j][n+1][k]=0.0;

          for(j=1;j<=SDIV-1;j++)
             for(n=1;n<=RNS_lmax;n++) {
                f_rho[j][n+1][k]=0.0;
                f_gama[j][n+1][k]=0.0;
             }
 
 
        n=0;
          for(j=2;j<=SDIV-1;j++)
             for(k=2;k<=SDIV-1;k++) {
                sk=s_gp[k];
                sj=s_gp[j];
                sk1=1.0-sk;
                sj1=1.0-sj;

                if(k<j) 
                  f_rho[j][n+1][k] = sj1/(sj*sk1*sk1);
                else     
                  f_rho[j][n+1][k] = 1.0/(sk*sk1);
             }
 
        j=SDIV;
          for(n=1;n<=RNS_lmax;n++)
             for(k=1;k<=SDIV;k++) {
                f_rho[j][n+1][k] = 0.0;
                f_gama[j][n+1][k] = 0.0;
             }

        k=SDIV;
          for(j=1;j<=SDIV;j++)
              for(n=1;n<=RNS_lmax;n++) {
                 f_rho[j][n+1][k] = 0.0;
                 f_gama[j][n+1][k] = 0.0;
              }
   }

  n=0;
   for(i=1;i<=MDIV;i++)
      P_2n[i][n+1]=legendre(2*n,mu[i]); 

   for(i=1;i<=MDIV;i++)
     for(n=1;n<=RNS_lmax;n++) {
      P_2n[i][n+1]    = legendre(2*n,mu[i]);
      P1_2n_1[i][n+1] = plgndr(2*n-1 ,1,mu[i]);
    }

  array_free(f2n,1,RNS_lmax+1,1,SDIV);       


  for(m=1;m<=MDIV;m++) { 
     sin_theta[m] = sqrt(1.0-mu[m]*mu[m]);  
     theta[m] = asin(sin_theta[m]);
     /* try alternative: theta[m] = acos(mu[m]); */
  }
 
  r_e = (*r_e_new);
 
  /* RESCALE Omega_e */

  (*Omega_e) *= r_e;


  /* START MAIN ITERATIVE LOOP */

  do {
      ITERATION_COUNT += 1;
      if ( ITERATION_COUNT > MAX_ITERATION ) {
          CCTK_VERROR("Max number of iteration %d reached ! Exiting.", MAX_ITERATION);
	  }

 
      /* Rescale potentials and construct arrays with the potentials along
       | the equatorial and polar directions.
      */        

      for(s=1;s<=SDIV;s++) {
         for(m=1;m<=MDIV;m++) {
            rho[s][m] /= SQ(r_e);
            gama[s][m] /= SQ(r_e); 
            alpha[s][m] /= SQ(r_e);
            omega[s][m] *= r_e;
         }
         rho_mu_0[s]=rho[s][1];     
         gama_mu_0[s]=gama[s][1];   
         omega_mu_0[s]=omega[s][1]; 
         rho_mu_1[s]=rho[s][MDIV];  
         gama_mu_1[s]=gama[s][MDIV];
      }


      /* Find location of maximum energy density */

      if( is_uniform_rotation) {

	if(energy[1][1]==0.0) {  /* AT FIRST ITERATION SET s_max AT CENTER */
          s_max=0.0;
        }else{
	      for(s=1;s<=SDIV/2-1;s++) {
                 d_e_ds_mu_0[s] = deriv_s(energy,s,1);
              }

              n_nearest = SDIV/4;
              s_max = interp(d_e_ds_mu_0, s_gp, SDIV/2-1, 0.0, &n_nearest);

	}

        if(s_max<0.0)
          s_max=0.0;
      }

 
      /* Compute new r_e. */ 

      r_e_old=r_e;
      r_p=r_ratio*r_e;                          
      s_p=r_p/(r_p+r_e);                        
  
      n_nearest= SDIV/2;
      gama_pole_h=interp(s_gp,gama_mu_1,SDIV,s_p,&n_nearest); 
      gama_equator_h=interp(s_gp,gama_mu_0,SDIV,s_e,&n_nearest);
      gama_center_h=gama[1][1];                    
  
      rho_pole_h=interp(s_gp,rho_mu_1,SDIV,s_p,&n_nearest);   
      rho_equator_h=interp(s_gp,rho_mu_0,SDIV,s_e,&n_nearest);
      rho_center_h=rho[1][1];                      
 
      r_e=sqrt(2*h_center/(gama_pole_h+rho_pole_h-gama_center_h-rho_center_h));


      /* Compute angular velocity Omega */
 
      if( is_uniform_rotation) {

        if(r_ratio==1.0) {
          Omega_h=0.0;
          omega_equator_h=0.0;
        } 
        else {
              omega_equator_h=interp(s_gp,omega_mu_0,SDIV,s_e, &n_nearest);
              term_in_Omega_h=1.0-exp(SQ(r_e)*(gama_pole_h+rho_pole_h
                                               -gama_equator_h-rho_equator_h));
              if(term_in_Omega_h>=0.0) 
                 Omega_h = omega_equator_h + exp(SQ(r_e)*rho_equator_h)
                                            *sqrt(term_in_Omega_h);
              else {
                    Omega_h=0.0;
                    printf("sqrt in Omega_h imaginary\n");
	      }
        }

	for(s=1;s<=2*SDIV/3;s++)
	  for(m=1;m<=MDIV;m++) 
	    Omega_diff[s][m] = Omega_h;

        /* this pacifies the compiler. The variable is not used in the uniform
         * case. */
        Omega_c = atof("nan");
      }
      else { 
              omega_equator_h=interp(s_gp,omega_mu_0,SDIV,s_e, &n_nearest);
 
	          if(n_of_it==0) {   
                     Omega_h = omega_equator_h + exp(SQ(r_e)*rho_equator_h)
                               *sqrt(1.0-exp(SQ(r_e)*(gama_pole_h+rho_pole_h
                               -gama_equator_h-rho_equator_h)));

                     Omega_e_lower_bracket = 0.0;
                     Omega_e_upper_bracket = 1.0*Omega_h;
                  }               
                  else {
                       Omega_e_lower_bracket = (*Omega_e)/2.0;

                       if(A_diff<(0.7/r_e)) {
                         Omega_e_upper_bracket = 2.0*(*Omega_e);
                       }
                       else {
			     if(A_diff<(0.8/r_e)) {
                               Omega_e_upper_bracket = 1.75*(*Omega_e);
                             } 
                             else {
                               Omega_e_upper_bracket = 1.5*(*Omega_e);
                             } 
		            }
                  }

		  tol = DBL_EPSILON;

                  (*Omega_e) = zbrent_diff( diff_rotation, r_e, rho_equator_h, 
                                   gama_equator_h, omega_equator_h, rho_pole_h,
                                   gama_pole_h, 
                                   Omega_e_lower_bracket, 
                                   Omega_e_upper_bracket, tol, A_diff);

  		  Omega_c = (*Omega_e) + SQ(1.0/A_diff)*( (*Omega_e)
                              -omega_equator_h)
                              *exp(-2.0*SQ(r_e)*rho_equator_h)/( 1.0 - 
                              SQ(( (*Omega_e) -omega_equator_h)*exp(-SQ(r_e)
                              *rho_equator_h)) );


                  Omega_h =  omega_equator_h + exp(SQ(r_e)*rho_equator_h)
                               *sqrt(1.0-exp(SQ(r_e)*(gama_pole_h+rho_pole_h
                               -gama_equator_h-rho_equator_h)));

                  s=1;
 	          for(m=1;m<=MDIV-1;m++) 
                          Omega_diff[s][m] = Omega_c;

                  m=MDIV;
                  for(s=1;s<=2*SDIV/3;s++) 
                          Omega_diff[s][m] = Omega_c;



                  for(s=2;s<=2*SDIV/3;s++) 
                     for(m=1;m<=MDIV-1;m++) { 

                       Omega_lower_bracket = 0.0;
                     
                       Omega_upper_bracket = 1.2*Omega_diff[s-1][m];
                       

                         tol = DBL_EPSILON;

                         Omega_diff[s][m] = zbrent_rot( rotation_law, r_e,
                                                       rho[s][m], omega[s][m],
                                                       s_gp[s],mu[m], Omega_c,
                                                       Omega_lower_bracket,
                                                       Omega_upper_bracket,
                                                       tol, A_diff);     
                     }
      }


      /* Compute velocity, energy density and pressure. */
 
      n_nearest=n_tab/2; 

      if( is_uniform_rotation) {
	
	/* uniform rotation */

        for(s=1;s<=SDIV;s++) {
           sgp=s_gp[s];

           for(m=1;m<=MDIV;m++) {
              rsm=rho[s][m];

              if(r_ratio==1.0 || s > (SDIV/2+2) ) 
                    velocity_sq[s][m]=0.0;
              else {
                    velocity_sq[s][m]=SQ((Omega_h-omega[s][m])*(sgp/(1.0-sgp))
                                      *sin_theta[m]*exp(-rsm*SQ(r_e)));
	      }
              if(velocity_sq[s][m]>=1.0) 
                velocity_sq[s][m]=0.0;
           }
	}

	/* differential rotation */

      } else {
            for(s=1;s<=SDIV;s++) {
	      sgp=s_gp[s];

	      for(m=1;m<=MDIV;m++) {
		rsm=rho[s][m];
		
		if(r_ratio==1.0 || s > (SDIV/2+2) ) 
		  velocity_sq[s][m]=0.0;
		else 
		  velocity_sq[s][m]=SQ( (Omega_diff[s][m]-omega[s][m])
					*(sgp/(1.0-sgp))*sin_theta[m]
					*exp(-rsm*SQ(r_e)) );
		
		if(velocity_sq[s][m]>=1.0) 
		  velocity_sq[s][m]=0.0;
	      }
	    }
      }
      
      
      for(s=1;s<=SDIV;s++) {
         sgp=s_gp[s];

         for(m=1;m<=MDIV;m++) {
            rsm=rho[s][m];
            
            if( is_uniform_rotation)
              enthalpy[s][m] = enthalpy_min + 0.5*(SQ(r_e)*(gama_pole_h
				+rho_pole_h-gama[s][m]-rsm) 
                               -log(1.0-velocity_sq[s][m]));
            else
              enthalpy[s][m] = enthalpy_min + 0.5*(SQ(r_e)*(gama_pole_h
                               +rho_pole_h-gama[s][m]-rsm)-log(1.0
                               -velocity_sq[s][m])+SQ(A_diff*(Omega_diff[s][m]
                               -Omega_c)) );
  
            if((enthalpy[s][m]<=enthalpy_min) || (sgp>s_e)) {
                  pressure[s][m]=0.0;
                  energy[s][m]=0.0; 
	    }
	    else { if(strcmp(eos_type,"tab")==0) {

	      pressure[s][m]=p_at_h(enthalpy[s][m], log_p_tab, 
				    log_h_tab, n_tab, &n_nearest);
	      energy[s][m]=e_at_p(pressure[s][m], log_e_tab, 
				  log_p_tab, n_tab, &n_nearest, eos_type, eos_k,
				  Gamma_P);

	    }
	    else { if(strcmp(eos_type,"poly")==0) {
		rho0sm=pow(((Gamma_P-1.0)/(Gamma_P*eos_k))
			 *(exp(enthalpy[s][m])-1.0),1.0/(Gamma_P-1.0));

	      pressure[s][m]=eos_k*pow(rho0sm,Gamma_P);

	      energy[s][m]=pressure[s][m]/(Gamma_P-1.0)+rho0sm;
	    }
	    else { /* uniform density - activate later */
	      /*
		pressure[s][m]=exp(enthalpy[s][m])-1.0;
		energy[s][m]=e_center;
	      */
	      printf(" NO NO NO NO NO OTHER EOS ");

	    }  
	    }


	    }  
	    
            /* Rescale back metric potentials (except omega) */
	    
            rho[s][m] *= SQ(r_e);
            gama[s][m] *= SQ(r_e);
            alpha[s][m] *= SQ(r_e);
	 }
      }


      /* Compute metric potentials */

      S_gama = array_allocate(1,SDIV,1,MDIV);
      S_rho = array_allocate(1,SDIV,1,MDIV);
      S_omega = array_allocate(1,SDIV,1,MDIV);


      if( is_uniform_rotation) {

       for(s=1;s<=SDIV;s++)
         for(m=1;m<=MDIV;m++) {
            rsm=rho[s][m];
            gsm=gama[s][m];
            omsm=omega[s][m];
            esm=energy[s][m];
            psm=pressure[s][m];
            e_gsm=exp(0.5*gsm);
            e_rsm=exp(-rsm);
            v2sm=velocity_sq[s][m];
            mum=mu[m];            
            m1=1.0-SQ(mum);
            sgp=s_gp[s];
            s_1=1.0-sgp;
            s1=sgp*s_1;
            s2=SQ(sgp/s_1);  

            ea=16.0*PI*exp(2.0*alpha[s][m])*SQ(r_e);
 
            if(s==1) {
              d_gama_s=0.0;
              d_gama_m=0.0;
              d_rho_s=0.0;
              d_rho_m=0.0;
              d_omega_s=0.0;
              d_omega_m=0.0;
            }else{
                 d_gama_s=deriv_s(gama,s,m);
                 d_gama_m=deriv_m(gama,s,m);
                 d_rho_s=deriv_s(rho,s,m);
                 d_rho_m=deriv_m(rho,s,m);
                 d_omega_s=deriv_s(omega,s,m);
                 d_omega_m=deriv_m(omega,s,m);
	     }

            S_rho[s][m] = e_gsm*(0.5*ea*(esm + psm)*s2*(1.0+v2sm)/(1.0-v2sm)
  
                          + s2*m1*SQ(e_rsm)*(SQ(s1*d_omega_s) 
                       
                          + m1*SQ(d_omega_m))
                         
                          + s1*d_gama_s - mum*d_gama_m + 0.5*rsm*(ea*psm*s2  
 
                          - s1*d_gama_s*(0.5*s1*d_gama_s+1.0) 
 
                          - d_gama_m*(0.5*m1*d_gama_m-mum)));

            S_gama[s][m] = e_gsm*(ea*psm*s2 + 0.5*gsm*(ea*psm*s2 - 0.5*SQ(s1

                           *d_gama_s) - 0.5*m1*SQ(d_gama_m)));

            S_omega[s][m]=e_gsm*e_rsm*( -ea*(Omega_h-omsm)*(esm+psm)

                          *s2/(1.0-v2sm) + omsm*( -0.5*ea*(((1.0+v2sm)*esm 
                           
                          + 2.0*v2sm*psm)/(1.0-v2sm))*s2 

                          - s1*(2*d_rho_s+0.5*d_gama_s)

                          + mum*(2*d_rho_m+0.5*d_gama_m) + 0.25*SQ(s1)*(4

                          *SQ(d_rho_s)-SQ(d_gama_s)) + 0.25*m1*(4*SQ(d_rho_m)

                          - SQ(d_gama_m)) - m1*SQ(e_rsm)*(SQ(SQ(sgp)*d_omega_s)

                          + s2*m1*SQ(d_omega_m))));
	 }
     }
     else{


       for(s=1;s<=SDIV;s++)
         for(m=1;m<=MDIV;m++) {
            rsm=rho[s][m];
            gsm=gama[s][m];
            omsm=omega[s][m];
            esm=energy[s][m];
            psm=pressure[s][m];
            e_gsm=exp(0.5*gsm);
            e_rsm=exp(-rsm);
            v2sm=velocity_sq[s][m];
            mum=mu[m];            
            m1=1.0-SQ(mum);
            sgp=s_gp[s];
            s_1=1.0-sgp;
            s1=sgp*s_1;
            s2=SQ(sgp/s_1);  

            ea=16.0*PI*exp(2.0*alpha[s][m])*SQ(r_e);
 
            if(s==1) {
              d_gama_s=0.0;
              d_gama_m=0.0;
              d_rho_s=0.0;
              d_rho_m=0.0;
              d_omega_s=0.0;
              d_omega_m=0.0;
            }else{
                 d_gama_s=deriv_s(gama,s,m);
                 d_gama_m=deriv_m(gama,s,m);
                 d_rho_s=deriv_s(rho,s,m);
                 d_rho_m=deriv_m(rho,s,m);
                 d_omega_s=deriv_s(omega,s,m);
                 d_omega_m=deriv_m(omega,s,m);
	     }

            S_rho[s][m] = e_gsm*(0.5*ea*(esm + psm)*s2*(1.0+v2sm)/(1.0-v2sm)
  
                          + s2*m1*SQ(e_rsm)*(SQ(s1*d_omega_s) + m1*SQ(d_omega_m))
                         
                          + s1*d_gama_s - mum*d_gama_m + 0.5*rsm*(ea*psm*s2  
 
                          - s1*d_gama_s*(0.5*s1*d_gama_s+1.0) 
 
                          - d_gama_m*(0.5*m1*d_gama_m-mum)));

            S_gama[s][m] = e_gsm*(ea*psm*s2 + 0.5*gsm*(ea*psm*s2 - 0.5*SQ(s1

                           *d_gama_s) - 0.5*m1*SQ(d_gama_m)));

            S_omega[s][m]=e_gsm*e_rsm*( -ea*(Omega_diff[s][m]-omsm)*(esm+psm)

                          *s2/(1.0-v2sm) + omsm*( -0.5*ea*(((1.0+v2sm)*esm + 2.0 
      
                          *v2sm*psm)/(1.0-v2sm))*s2 - s1*(2*d_rho_s+0.5*d_gama_s)

                          + mum*(2*d_rho_m+0.5*d_gama_m) + 0.25*SQ(s1)*(4

                          *SQ(d_rho_s)-SQ(d_gama_s)) + 0.25*m1*(4*SQ(d_rho_m)

                          - SQ(d_gama_m)) - m1*SQ(e_rsm)*(SQ(SQ(sgp)*d_omega_s) 

                          + s2*m1*SQ(d_omega_m))));
	 }
     }
    


      /* ANGULAR INTEGRATION */
   
      D1_rho = array_allocate(1,RNS_lmax+1,1,SDIV);  
      D1_gama = array_allocate(1,RNS_lmax+1,1,SDIV); 
      D1_omega = array_allocate(1,RNS_lmax+1,1,SDIV);

      n=0;
      for(k=1;k<=SDIV;k++) {      

         for(m=1;m<=MDIV-2;m+=2) {
               sum_rho += (DM/3.0)*(P_2n[m][n+1]*S_rho[k][m]
                          + 4.0*P_2n[m+1][n+1]*S_rho[k][m+1] 
                          + P_2n[m+2][n+1]*S_rho[k][m+2]);
	 }

         D1_rho[n+1][k]=sum_rho;
         D1_gama[n+1][k]=0.0;
         D1_omega[n+1][k]=0.0;
         sum_rho=0.0;

      }

      for(n=1;n<=RNS_lmax;n++)
         for(k=1;k<=SDIV;k++) {      
            for(m=1;m<=MDIV-2;m+=2) {

               sum_rho += (DM/3.0)*(P_2n[m][n+1]*S_rho[k][m]
                          + 4.0*P_2n[m+1][n+1]*S_rho[k][m+1] 
                          + P_2n[m+2][n+1]*S_rho[k][m+2]);
                       
               sum_gama += (DM/3.0)*(sin((2.0*n-1.0)*theta[m])*S_gama[k][m]
                           +4.0*sin((2.0*n-1.0)*theta[m+1])*S_gama[k][m+1]
                           +sin((2.0*n-1.0)*theta[m+2])*S_gama[k][m+2]);
  
               sum_omega += (DM/3.0)*(sin_theta[m]*P1_2n_1[m][n+1]*S_omega[k][m]
                            +4.0*sin_theta[m+1]*P1_2n_1[m+1][n+1]*S_omega[k][m+1]
                            +sin_theta[m+2]*P1_2n_1[m+2][n+1]*S_omega[k][m+2]);
	    }
            D1_rho[n+1][k]=sum_rho;
            D1_gama[n+1][k]=sum_gama;
            D1_omega[n+1][k]=sum_omega;
            sum_rho=0.0;
            sum_gama=0.0;
            sum_omega=0.0;
	}


      array_free(S_gama,1,SDIV,1,MDIV);   
      array_free(S_rho,1,SDIV,1,MDIV);    
      array_free(S_omega,1,SDIV,1,MDIV);  



      /* RADIAL INTEGRATION */



      D2_rho = array_allocate(1,SDIV,1,RNS_lmax+1);  
      D2_gama = array_allocate(1,SDIV,1,RNS_lmax+1); 
      D2_omega = array_allocate(1,SDIV,1,RNS_lmax+1);


      n=0;
      for(s=1;s<=SDIV;s++) {
            for(k=1;k<=SDIV-2;k+=2) { 
               sum_rho += (DS/3.0)*( f_rho[s][n+1][k]*D1_rho[n+1][k] 
                          + 4.0*f_rho[s][n+1][k+1]*D1_rho[n+1][k+1]
                          + f_rho[s][n+1][k+2]*D1_rho[n+1][k+2]);
 	    }
	    D2_rho[s][n+1]=sum_rho;
	    D2_gama[s][n+1]=0.0;
	    D2_omega[s][n+1]=0.0;
            sum_rho=0.0;
	 }

 
      for(s=1;s<=SDIV;s++)
         for(n=1;n<=RNS_lmax;n++) {
            for(k=1;k<=SDIV-2;k+=2) { 
               sum_rho += (DS/3.0)*( f_rho[s][n+1][k]*D1_rho[n+1][k] 
                          + 4.0*f_rho[s][n+1][k+1]*D1_rho[n+1][k+1]
                          + f_rho[s][n+1][k+2]*D1_rho[n+1][k+2]);
 
               sum_gama += (DS/3.0)*( f_gama[s][n+1][k]*D1_gama[n+1][k] 
                           + 4.0*f_gama[s][n+1][k+1]*D1_gama[n+1][k+1]
                           + f_gama[s][n+1][k+2]*D1_gama[n+1][k+2]);
     
               if(k<s && k+2<=s) 
                 sum_omega += (DS/3.0)*( f_rho[s][n+1][k]*D1_omega[n+1][k] 
                              + 4.0*f_rho[s][n+1][k+1]*D1_omega[n+1][k+1]
                              + f_rho[s][n+1][k+2]*D1_omega[n+1][k+2]);
               else {
                 if(k>=s) 
                   sum_omega += (DS/3.0)*( f_gama[s][n+1][k]*D1_omega[n+1][k] 
                                + 4.0*f_gama[s][n+1][k+1]*D1_omega[n+1][k+1]
                                + f_gama[s][n+1][k+2]*D1_omega[n+1][k+2]);
                 else
                   sum_omega += (DS/3.0)*( f_rho[s][n+1][k]*D1_omega[n+1][k] 
                                + 4.0*f_rho[s][n+1][k+1]*D1_omega[n+1][k+1]
                                + f_gama[s][n+1][k+2]*D1_omega[n+1][k+2]);
               }
	    }
	    D2_rho[s][n+1]=sum_rho;
	    D2_gama[s][n+1]=sum_gama;
	    D2_omega[s][n+1]=sum_omega;
            sum_rho=0.0;
            sum_gama=0.0;
            sum_omega=0.0;
	 }
 
      array_free(D1_rho,1,RNS_lmax+1,1,SDIV);   
      array_free(D1_gama,1,RNS_lmax+1,1,SDIV);  
      array_free(D1_omega,1,RNS_lmax+1,1,SDIV); 


      /* SUMMATION OF COEFFICIENTS */

      for(s=1;s<=SDIV;s++) 
         for(m=1;m<=MDIV;m++) {

            gsm=gama[s][m];
            rsm=rho[s][m];
            omsm=omega[s][m];             
            e_gsm=exp(-0.5*gsm);
            e_rsm=exp(rsm);
            temp1=sin_theta[m];

            sum_rho += -e_gsm*P_2n[m][0+1]*D2_rho[s][0+1]; 

            for(n=1;n<=RNS_lmax;n++) {

               sum_rho += -e_gsm*P_2n[m][n+1]*D2_rho[s][n+1]; 

               if(m==MDIV) {             
                 sum_omega += 0.5*e_rsm*e_gsm*D2_omega[s][n+1]; 
                 sum_gama += -(2.0/PI)*e_gsm*D2_gama[s][n+1];   
	       }
               else { 
                     sum_omega += -e_rsm*e_gsm*(P1_2n_1[m][n+1]/(2.0*n
                                  *(2.0*n-1.0)*temp1))*D2_omega[s][n+1];
  
                     sum_gama += -(2.0/PI)*e_gsm*(sin((2.0*n-1.0)*theta[m])
                                 /((2.0*n-1.0)*temp1))*D2_gama[s][n+1];   
	       }
	    }
	   
            rho[s][m]=rsm + cf*(sum_rho-rsm);
            gama[s][m]=gsm + cf*(sum_gama-gsm);
            omega[s][m]=omsm + cf*(sum_omega-omsm);

            sum_omega=0.0;
            sum_rho=0.0;
            sum_gama=0.0; 
	  }

      array_free(D2_rho,1,SDIV,1,RNS_lmax+1);   
      array_free(D2_gama,1,SDIV,1,RNS_lmax+1);  
      array_free(D2_omega,1,SDIV,1,RNS_lmax+1); 


      /* CHECK FOR DIVERGENCE */

      if(fabs(omega[2][1])>100.0 || fabs(rho[2][1])>100.0 
         || fabs(gama[2][1])>300.0) {
         a_check=200; 
         break;
      }


      /* TREAT SPHERICAL CASE */

      if(r_ratio==1.0) {
        for(s=1;s<=SDIV;s++)
           for(m=1;m<=MDIV;m++) {
              rho[s][m]=rho[s][1];
              gama[s][m]=gama[s][1];
              omega[s][m]=0.0;          
	   }
      }


      /* TREAT INFINITY WHEN SMAX=1.0 */

      if(SMAX==1.0) {
         for(m=1;m<=MDIV;m++) {
            rho[SDIV][m]=0.0;
            gama[SDIV][m]=0.0;
            omega[SDIV][m]=0.0;
	 }
      } 


      /* COMPUTE FIRST ORDER DERIVATIVES OF GAMA */


      da_dm = array_allocate(1,SDIV,1,MDIV);
      dgds = array_allocate(1,SDIV,1,MDIV);
      dgdm = array_allocate(1,SDIV,1,MDIV);

      for(s=1;s<=SDIV;s++)
         for(m=1;m<=MDIV;m++) {
            dgds[s][m]=deriv_s(gama,s,m);
            dgdm[s][m]=deriv_m(gama,s,m);
	 }



      /* ALPHA */
 
      if(r_ratio==1.0) {
        for(s=1;s<=SDIV;s++)
           for(m=1;m<=MDIV;m++)
              da_dm[s][m]=0.0;
      } 
      else {
            for(s=2;s<=SDIV;s++)
               for(m=1;m<=MDIV;m++) {

                  da_dm[1][m]=0.0; 
       
                  sgp=s_gp[s];
                  s1=sgp*(1.0-sgp);
                  mum=mu[m]; 
                  m1=1.0-SQ(mum);
          
                  d_gama_s=dgds[s][m];
                  d_gama_m=dgdm[s][m];
                  d_rho_s=deriv_s(rho,s,m);
                  d_rho_m=deriv_m(rho,s,m);
                  d_omega_s=deriv_s(omega,s,m);
                  d_omega_m=deriv_m(omega,s,m);
                  d_gama_ss=s1*deriv_s(dgds,s,m)+(1.0-2.0*sgp)
                                               *d_gama_s;
                  d_gama_mm=m1*deriv_m(dgdm,s,m)-2.0*mum*d_gama_m;  
                  d_gama_sm=deriv_sm(gama,s,m);

           temp1=2.0*SQ(sgp)*(sgp/(1.0-sgp))*m1*d_omega_s*d_omega_m

                *(1.0+s1*d_gama_s) - (SQ(SQ(sgp)*d_omega_s) - 
 
                SQ(sgp*d_omega_m/(1.0-sgp))*m1)*(-mum+m1*d_gama_m); 
  
           temp2=1.0/(m1 *SQ(1.0+s1*d_gama_s) + SQ(-mum+m1*d_gama_m));

           temp3=s1*d_gama_ss + SQ(s1*d_gama_s);
  
           temp4=d_gama_m*(-mum+m1*d_gama_m);
   
           temp5=(SQ(s1*(d_rho_s+d_gama_s)) - m1*SQ(d_rho_m+d_gama_m))

                 *(-mum+m1*d_gama_m);

           temp6=s1*m1*(0.5*(d_rho_s+d_gama_s)* (d_rho_m+d_gama_m) 
  
                + d_gama_sm + d_gama_s*d_gama_m)*(1.0 + s1*d_gama_s); 

           temp7=s1*mum*d_gama_s*(1.0+s1*d_gama_s);

           temp8=m1*exp(-2*rho[s][m]);

          da_dm[s][m] = -0.5*(d_rho_m+d_gama_m) - temp2*(0.5*(temp3 - 

            d_gama_mm - temp4)*(-mum+m1*d_gama_m) + 0.25*temp5 

            - temp6 +temp7 + 0.25*temp8*temp1);	 
       }
   }

      for(s=1;s<=SDIV;s++) {
         alpha[s][1]=0.0;
         for(m=1;m<=MDIV-1;m++) 
            alpha[s][m+1]=alpha[s][m]+0.5*DM*(da_dm[s][m+1]+
                          da_dm[s][m]);
      }


   array_free(da_dm,1,SDIV,1,MDIV);
   array_free(dgds,1,SDIV,1,MDIV);
   array_free(dgdm,1,SDIV,1,MDIV);


      for(s=1;s<=SDIV;s++)
         for(m=1;m<=MDIV;m++) {
            alpha[s][m] += -alpha[s][MDIV]+0.5*(gama[s][MDIV]-rho[s][MDIV]);

            if(alpha[s][m]>=300.0) {
              a_check=200;
              break;
            }
            omega[s][m] /= r_e;
         }

      if(SMAX==1.0) {
         for(m=1;m<=MDIV;m++)
            alpha[SDIV][m] = 0.0;
      }

      if(a_check==200)
        break;

      dif=fabs(r_e_old-r_e)/r_e;
      n_of_it++;

      if(print_dif!=0)
         printf(" %4.3e\n",dif);

  } while(dif> accuracy || n_of_it<=2);   /* end do */



  /* COMPUTE OMEGA */

  if(strcmp(eos_type,"tab")==0)
/*      (*Omega) = Omega_h*C/(r_e*sqrt(KAPPA)); */ /*this should be a conversion*/
      (*Omega) = Omega_h/r_e;
  else
    (*Omega) = Omega_h/r_e;


    /* Rescale back Omega_diff[s][m] */

     /*
    if( strcmp(rotation_type,"diff")==0) {
     */
         for(s=1; s<=2*SDIV/3; s++)
	    for(m=1; m<=MDIV; m++) {
               Omega_diff[s][m] /= r_e;
            }
	 /*
    }
	 */

    /* RESCALE BACK Omega_e */

    (*Omega_e) /= r_e;


    /* UPDATE r_e_new */

    (*r_e_new) = r_e;


    tensor_free(f_rho, 1,SDIV,1,RNS_lmax+1,1,SDIV);
    tensor_free(f_gama,1,SDIV,1,RNS_lmax+1,1,SDIV);
    array_free(P_2n,   1,MDIV,1,RNS_lmax+1);
    array_free(P1_2n_1,1,MDIV,1,RNS_lmax+1);
}
