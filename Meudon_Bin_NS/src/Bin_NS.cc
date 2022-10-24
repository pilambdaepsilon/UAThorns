/* (c) 2009 Erik Schnetter
 * (c) 2010 Frank Loeffler */

#include <cstdio>
#include <cstring>
#include <cassert>
#include <vector>
#include <ios>
#include <iostream>
#include <stdlib.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <bin_ns.h>
#include <unites.h>
#define MAX_NTAB 16001
#define IMAX(a,b) ( a>b ? a : b ) 
#define IMIN(a,b) ( a<b ? a : b ) 
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-15
#endif

/*C*/
/***************************************************************************/
/* Routine that locates nearest grid point for a given value.              */
/***************************************************************************/
void hunt(double xx[], int n, double x, int *jlo)
{ 
	int jm,jhi,inc,ascnd;

	ascnd=(xx[n] > xx[1]);
	if (*jlo <= 0 || *jlo > n) {
		*jlo=0;
		jhi=n+1;
	} else {
		inc=1;
		if (x >= xx[*jlo] == ascnd) {
			if (*jlo == n) return;
			jhi=(*jlo)+1;
			while (x >= xx[jhi] == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n) {
					jhi=n+1;
					break;
				}
			}
		} else {
			if (*jlo == 1) {
				*jlo=0;
				return;
			}
			jhi=(*jlo);
			*jlo -= 1;
			while (x < xx[*jlo] == ascnd) {
				jhi=(*jlo);
				inc += inc;
				*jlo=jhi-inc;
				if (*jlo < 1) {
					*jlo=0;
					break;
				}
			}
		}
	}
	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if (x > xx[jm] == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
}
/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate using four points around xb.       */  
/*************************************************************************/
double interp(double xp[], 
              double yp[], 
              int    np ,
              double xb, 
              int    *n_nearest_pt)
{ 
 int k,        /* index of 1st point */
     m=4;      /* degree of interpolation */ 
 
 double y;     /* intermediate value */

 hunt(xp,np,xb,n_nearest_pt);

 k=IMIN(IMAX((*n_nearest_pt)-(m-1)/2,1),np+1-m);

 if( xb==xp[k] ||  xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) 
    xb += DBL_EPSILON;

 y= (xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
        ((xp[k]-xp[k+1])*(xp[k]-xp[k+2])*(xp[k]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/
       ((xp[k+1]-xp[k])*(xp[k+1]-xp[k+2])*(xp[k+1]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/
       ((xp[k+2]-xp[k])*(xp[k+2]-xp[k+1])*(xp[k+2]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
       ((xp[k+3]-xp[k])*(xp[k+3]-xp[k+1])*(xp[k+3]-xp[k+2]));

 return (y);
}
/*************************************************************************/
/* Load Beta equil file.                                                        */
/*************************************************************************/
void load_beta_equil( const char beta_equil_file[],
               double log_rho0_table[MAX_NTAB],
               double Y_e_table[MAX_NTAB],
               int *n_tab_beta)
{
 int i;                    /* counter */

 /*constants to convert from cgs to cactus units c=G=M_sun=1.0 */
 CCTK_REAL const cactusM= (5.028916268544129e-34);    /*  1/g  */
 CCTK_REAL const cactusL= (6.772400341316594e-06);    /*  1/cm */
 CCTK_REAL const cactusT= (2.0303145448833407e5);    /*  1/s  */
 CCTK_REAL const cactusV= (1.0/(cactusL*cactusL*cactusL));
 double rho0,               /* density */
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
       fscanf(f_beta,"%lf %lf \n",&rho0,&ye) ;
       log_rho0_table[i]=log10(rho0*cactusM*cactusV);     /* multiply by C^2 to get energy density */
       if(ye <= 0.036){
           Y_e_table[i] = 0.036;
       }
       else{
           Y_e_table[i]=ye;
       }
    }
}

using namespace std;

// define namespace here for old versions of Lorene that don't do so
namespace Lorene {}
using namespace Lorene;

int n_tab_beta;
double Y_e_tab[MAX_NTAB], log_rho0_tab_beta[MAX_NTAB];
int n_nearest_beta;
static void set_dt_from_domega (CCTK_ARGUMENTS,
                                CCTK_REAL const* const var,
                                CCTK_REAL      * const dtvar,
                                CCTK_REAL const& omega)
{
  DECLARE_CCTK_ARGUMENTS;

  int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  vector<CCTK_REAL> dxvar(npoints), dyvar(npoints);

  Diff_gv (cctkGH, 0, var, &dxvar[0], -1);
  Diff_gv (cctkGH, 1, var, &dyvar[0], -1);

#pragma omp parallel for
  for (int i=0; i<npoints; ++i) {
    CCTK_REAL const ephix = +y[i];
    CCTK_REAL const ephiy = -x[i];
    CCTK_REAL const dphi_var = ephix * dxvar[i] + ephiy * dyvar[i];
    dtvar[i] = omega * dphi_var;
  }
}

extern "C"
void Meudon_Bin_NS_initialise (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO ("Setting up LORENE Bin_NS initial data");
  if(init_real){CCTK_INFO ("(with realistic EOS)"); load_beta_equil(beta_file, log_rho0_tab_beta, Y_e_tab, &n_tab_beta);}


  // Meudon data are distributed in SI units (MKSA).  Here are some
  // conversion factors.
  // Be aware: these are the constants Lorene uses. They do differ from other
  // conventions, but they gave the best results in some tests.

  CCTK_REAL const c_light  = Unites::c_si;      // speed of light [m/s]
  CCTK_REAL const nuc_dens = Unites::rhonuc_si; // Nuclear density as used in Lorene units [kg/m^3]
  CCTK_REAL const G_grav   = Unites::g_si;      // gravitational constant [m^3/kg/s^2]
  CCTK_REAL const M_sun    = Unites::msol_si;   // solar mass [kg]

  // Cactus units in terms of SI units:
  // (These are derived from M = M_sun, c = G = 1, and using 1/M_sun
  // for the magnetic field)
  CCTK_REAL const cactusM = M_sun;
  CCTK_REAL const cactusL = cactusM * G_grav / pow(c_light,2);
  CCTK_REAL const cactusT = cactusL / c_light;

  // Other quantities in terms of Cactus units
  CCTK_REAL const coord_unit = cactusL / 1.0e+3;         // from km (~1.477)
  CCTK_REAL const rho_unit   = cactusM / pow(cactusL,3); // from kg/m^3
  CCTK_INT keyerr = 0, anyerr = 0;

  CCTK_INFO ("Setting up coordinates");

  int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  vector<double> xx(npoints), yy(npoints), zz(npoints);

  //Get EOS_Omni handle 

  if (!(*init_eos_key = EOS_Omni_GetHandle(eos_table)))
    CCTK_WARN(0,"Cannot get initial eos handle, aborting...");
  
  CCTK_VInfo(CCTK_THORNSTRING, "Meudon_Bin_NS will use the %s equation of state.", eos_table);
  CCTK_VInfo(CCTK_THORNSTRING, "Meudon_Bin_NS will use the %d eos handle", *init_eos_key);

 


#pragma omp parallel for
  for (int i=0; i<npoints; ++i) {
    xx[i] = x[i] * coord_unit;
    yy[i] = y[i] * coord_unit;
    zz[i] = z[i] * coord_unit;
  }

  // --------------------------------------------------------------
  //   CHECKING FILE NAME EXISTENCE
  // --------------------------------------------------------------
  FILE *file;
  if ((file = fopen(filename, "r")) != NULL) 
     fclose(file);
  else {
     CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                 "File \"%s\" does not exist. ABORTING", filename);
  }
  // Handle potentially different EOS table directory. LORENE recieves that via
  // environment variable
  if (strlen(eos_table_filepath) > 0) {
    if (setenv("LORENE_TABULATED_EOS_PATH", eos_table_filepath, 1)) {
      CCTK_ERROR("Unable to set environment variable LORENE_TABULATED_EOS_PATH");

    }
  }


  CCTK_VInfo (CCTK_THORNSTRING, "Reading from file \"%s\"", filename);

  try {
  Bin_NS bin_ns (npoints, &xx[0], &yy[0], &zz[0], filename);

  CCTK_VInfo (CCTK_THORNSTRING, "omega [rad/s]:       %g", bin_ns.omega);
  CCTK_VInfo (CCTK_THORNSTRING, "dist [km]:           %g", bin_ns.dist);
  CCTK_VInfo (CCTK_THORNSTRING, "dist_mass [km]:      %g", bin_ns.dist_mass);
  CCTK_VInfo (CCTK_THORNSTRING, "mass1_b [M_sun]:     %g", bin_ns.mass1_b);
  CCTK_VInfo (CCTK_THORNSTRING, "mass2_b [M_sun]:     %g", bin_ns.mass2_b);
  CCTK_VInfo (CCTK_THORNSTRING, "mass_ADM [M_sun]:    %g", bin_ns.mass_adm);
  CCTK_VInfo (CCTK_THORNSTRING, "L_tot [G M_sun^2/c]: %g", bin_ns.angu_mom);
  CCTK_VInfo (CCTK_THORNSTRING, "rad1_x_comp [km]:    %g", bin_ns.rad1_x_comp);
  CCTK_VInfo (CCTK_THORNSTRING, "rad1_y [km]:         %g", bin_ns.rad1_y);
  CCTK_VInfo (CCTK_THORNSTRING, "rad1_z [km]:         %g", bin_ns.rad1_z);
  CCTK_VInfo (CCTK_THORNSTRING, "rad1_x_opp [km]:     %g", bin_ns.rad1_x_opp);
  CCTK_VInfo (CCTK_THORNSTRING, "rad2_x_comp [km]:    %g", bin_ns.rad2_x_comp);
  CCTK_VInfo (CCTK_THORNSTRING, "rad2_y [km]:         %g", bin_ns.rad2_y);
  CCTK_VInfo (CCTK_THORNSTRING, "rad2_z [km]:         %g", bin_ns.rad2_z);
  CCTK_VInfo (CCTK_THORNSTRING, "rad2_x_opp [km]:     %g", bin_ns.rad2_x_opp);
  // LORENE's EOS is in terms on number density n = rho/m_nucleon:
  // P = K n^Gamma
  // to convert to SI units:
  // K_SI(n) = K_LORENE rho_nuc c^2 / n_nuc^gamma
  // Converting this to be in terms of the mass density rho = n m_nucleon gets
  // changes n_nuc to rho_nuc:
  // K_SI(rho) = K_LORENE c^2 / rho_nuc^(gamma-1)
  // In SI units P has units of M / (L T^2) and rho has units of M/L^3 thus
  // K_SI has units of (L^3/M)^Gamma M/(L T^2).
  // In Cactus units P and rho have the same units thus K_Cactus is unitless.
  // Conversion between K_SI and K_Cactus thus amounts to dividing out the
  // units of the SI quantity.
  double K = bin_ns.kappa_poly1 * pow((pow(c_light, 6.0) /
             ( pow(G_grav, 3.0) * M_sun * M_sun *
               nuc_dens )),bin_ns.gamma_poly1-1.);
  CCTK_VInfo (CCTK_THORNSTRING, "K [ET unit]:         %.15g", K);

  assert (bin_ns.np == npoints);

  CCTK_INFO ("Filling in Cactus grid points");

#pragma omp parallel for
  for (int i=0; i<npoints; ++i) {

    if (CCTK_EQUALS(initial_lapse, "Meudon_Bin_NS")) { 
      alp[i] = bin_ns.nnn[i];
    }

    if (CCTK_EQUALS(initial_shift, "Meudon_Bin_NS")) { 
      betax[i] = -bin_ns.beta_x[i];
      betay[i] = -bin_ns.beta_y[i];
      betaz[i] = -bin_ns.beta_z[i];
    }

    if (CCTK_EQUALS(initial_data, "Meudon_Bin_NS")) {
      gxx[i] = bin_ns.g_xx[i];
      gxy[i] = bin_ns.g_xy[i];
      gxz[i] = bin_ns.g_xz[i];
      gyy[i] = bin_ns.g_yy[i];
      gyz[i] = bin_ns.g_yz[i];
      gzz[i] = bin_ns.g_zz[i];

      kxx[i] = bin_ns.k_xx[i] * coord_unit;
      kxy[i] = bin_ns.k_xy[i] * coord_unit;
      kxz[i] = bin_ns.k_xz[i] * coord_unit;
      kyy[i] = bin_ns.k_yy[i] * coord_unit;
      kyz[i] = bin_ns.k_yz[i] * coord_unit;
      kzz[i] = bin_ns.k_zz[i] * coord_unit;
    }

    if (CCTK_EQUALS(initial_data, "Meudon_Bin_NS")) {
      rho[i] = bin_ns.nbar[i] / rho_unit;
      if(init_real){
      if(rho[i]>=1e-7){
          double yeres = interp(log_rho0_tab_beta, Y_e_tab, n_tab_beta,log10(rho[i]), &n_nearest_beta);
	  if(yeres <=0.036){yeres = 0.036;} 
	  Y_e[i] = yeres;
          temperature[i] = 0.1;
      }
      else{Y_e[i] = 0.25;}
      temperature[i] = 0.1;
      }
      if (!recalculate_eps){ //we don't know the temperature, so assume epsilon from ID is correct.
        eps[i] = bin_ns.ener_spec[i];
      }
      // Pressure from EOS_Omni call 
      if (CCTK_ActiveTimeLevelsVN(cctkGH, "HydroBase::temperature") > 0 &&
          CCTK_ActiveTimeLevelsVN(cctkGH, "HydroBase::Y_e") > 0)
      {
        EOS_Omni_press(*init_eos_key,recalculate_eps,eos_precision,1,&(rho[i]),&(eps[i]),
                       &(temperature[i]),&(Y_e[i]),&(press[i]),&keyerr,&anyerr);


      }
      else
      {
        EOS_Omni_press(*init_eos_key,recalculate_eps,eos_precision,1,&(rho[i]),&(eps[i]),
                       NULL,NULL,&(press[i]),&keyerr,&anyerr);
      }

      vel[i          ] = bin_ns.u_euler_x[i];
      vel[i+  npoints] = bin_ns.u_euler_y[i];
      vel[i+2*npoints] = bin_ns.u_euler_z[i];

      // Especially the velocity is set to strange values outside of the
      // matter region, so take care of this in the following way
      if (rho[i] < 1.e-15) {
        rho[i          ] = 1.e-15;
        vel[i          ] = 0.0;
        vel[i+  npoints] = 0.0;
        vel[i+2*npoints] = 0.0;
        eps[i          ] = K * pow(rho[i], bin_ns.gamma_poly1-1.) / (bin_ns.gamma_poly1-1.);
        press[i        ] = K * pow(rho[i], bin_ns.gamma_poly1);
      }
    }

  } // for i

  {
    // Angular velocity
    CCTK_REAL const omega = bin_ns.omega * cactusT;

    // These initial data assume a helical Killing vector field

    if (CCTK_EQUALS(initial_lapse, "Meudon_Bin_NS")) { 
      if (CCTK_EQUALS (initial_dtlapse, "Meudon_Bin_NS")) {
        CCTK_INFO ("Calculating time derivatives of lapse");
        set_dt_from_domega (CCTK_PASS_CTOC, alp, dtalp, omega);
      } else if (CCTK_EQUALS (initial_dtlapse, "none") or CCTK_EQUALS(initial_dtlapse,"zero")) {
        // do nothing
      } else {
        CCTK_WARN (CCTK_WARN_ABORT, "internal error");
      }
    }

    if (CCTK_EQUALS(initial_shift, "Meudon_Bin_NS")) { 
      if (CCTK_EQUALS (initial_dtshift, "Meudon_Bin_NS")) {
        CCTK_INFO ("Calculating time derivatives of shift");
        set_dt_from_domega (CCTK_PASS_CTOC, betax, dtbetax, omega);
        set_dt_from_domega (CCTK_PASS_CTOC, betay, dtbetay, omega);
        set_dt_from_domega (CCTK_PASS_CTOC, betaz, dtbetaz, omega);
      } else if (CCTK_EQUALS (initial_dtshift, "none") or CCTK_EQUALS(initial_dtshift,"zero")) {
        // do nothing
      } else {
        CCTK_WARN (CCTK_WARN_ABORT, "internal error");
      }
    }
  }

  CCTK_INFO ("Done.");
  } catch (ios::failure e) {
    CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not read initial data from file '%s': %s", filename, e.what());
  }
}

