/* Changed by N. Stergioulas, 24/10/2001: include latest source files from cactus 3 */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "SpaceMask.h"

#include "cctk_DefineThorn.h"
#include "rnsid.h"

#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velx_p (&vel_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p (&vel_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p (&vel_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velx_p_p (&vel_p_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p_p (&vel_p_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p_p (&vel_p_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

#define sx (&scon[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sy (&scon[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sz (&scon[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sx_p (&scon_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sy_p (&scon_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sz_p (&scon_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sx_p_p (&scon_p_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sy_p_p (&scon_p_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sz_p_p (&scon_p_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])


static const char *rcsid="$Header$";
CCTK_FILEVERSION(Hydro_RNSID_rnsid_rfr_c)

void rnsid_init(CCTK_ARGUMENTS);

void rnsid_init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL  Omega,                    /* angular velocity */
    R_e,                               /* circumferential radius */
    r_e,                               /* coordinate radius */
    mass0;                             /* rest mass */

  CCTK_REAL eos_k, eos_ideal_fluid_gamma, rnsid_rho_min;

  CCTK_REAL *x_coord=0, *y_coord=0, *z_coord=0;

  CCTK_REAL gamma_center;  /* central value of 3-determinant of metric,
                              needed for applying perturbations */

  if ( (RNS_K < 0.0) || (RNS_Gamma < 0.0)) {
    CCTK_WARN(0,"RNS_K and RNS_Gamma must be greater than 0: using 100.0 and 2!");
    eos_k                 = 100.0;
    eos_ideal_fluid_gamma = 2.0;
  } else {
    eos_k                 = RNS_K;
    eos_ideal_fluid_gamma = RNS_Gamma;
  }

  rnsid_rho_min = RNS_rho_min;


 /* SET-UP INITIAL DATA */

  x_coord = x;
  y_coord = y;
  z_coord = z;

  hydro_rnsid(cctkGH, x_coord, y_coord, z_coord,eos_k, eos_ideal_fluid_gamma, rnsid_rho_min,
              &Omega, &R_e, &r_e, &mass0, &gamma_center);


  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::metric") > 1)
    {
      #pragma omp parallel for
      for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	{
	  gxx_p[i] = gxx[i];
	  gyy_p[i] = gyy[i];
	  gzz_p[i] = gzz[i];
	  gxy_p[i] = gxy[i];
	  gxz_p[i] = gxz[i];
	  gyz_p[i] = gyz[i];
	}
      if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::metric") > 2)
	{
	  for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	    {
	      gxx_p_p[i] = gxx[i];
	      gyy_p_p[i] = gyy[i];
	      gzz_p_p[i] = gzz[i];
	      gxy_p_p[i] = gxy[i];
	      gxz_p_p[i] = gxz[i];
	      gyz_p_p[i] = gyz[i];
	    }
	}
    }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::curv") > 1)
    {
      #pragma omp parallel for
      for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	{
	  kxx_p[i] = kxx[i];
	  kyy_p[i] = kyy[i];
	  kzz_p[i] = kzz[i];
	  kxy_p[i] = kxy[i];
	  kxz_p[i] = kxz[i];
	  kyz_p[i] = kyz[i];
	}
      if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::curv") > 2)
	{
          #pragma omp parallel for
	  for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	    {
	      kxx_p_p[i] = kxx[i];
	      kyy_p_p[i] = kyy[i];
	      kzz_p_p[i] = kzz[i];
	      kxy_p_p[i] = kxy[i];
	      kxz_p_p[i] = kxz[i];
	      kyz_p_p[i] = kyz[i];
	    }
	}
    }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::lapse") > 1)
    {
      #pragma omp parallel for
      for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	{
	  alp_p[i] = alp[i];
	}
      if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::lapse") > 2)
	{
          #pragma omp parallel for
	  for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	    {
	      alp_p_p[i] = alp[i];
	    }
	}
    }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::shift") > 1)
    {
      #pragma omp parallel for
      for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	{
	  betax_p[i] = betax[i];
	  betay_p[i] = betay[i];
	  betaz_p[i] = betaz[i];
	}
      if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::shift") > 2)
	{
          #pragma omp parallel for
	  for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	    {
	      betax_p_p[i] = betax[i];
	      betay_p_p[i] = betay[i];
	      betaz_p_p[i] = betaz[i];
	    }
	}
    }

  if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::rho") > 1)
    {
      #pragma omp parallel for
      for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	{
	  rho_p[i] = rho[i];
	  velx_p[i] = velx[i];
	  vely_p[i] = vely[i];
	  velz_p[i] = velz[i];
	  press_p[i] = press[i];
	  eps_p[i] = eps[i];
	  w_lorentz_p[i] = w_lorentz[i];
	}
      if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::rho") > 2)
	{
          #pragma omp parallel for
	  for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	    {
	      rho_p_p[i] = rho[i];
	      velx_p_p[i] = velx[i];
	      vely_p_p[i] = vely[i];
	      velz_p_p[i] = velz[i];
	      press_p_p[i] = press[i];
	      eps_p_p[i] = eps[i];
	      w_lorentz_p_p[i] = w_lorentz[i];
	    }
	}
    }

  return;
}

void RNSID_CheckParameters(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (timelevels < 2)
  {
      CCTK_PARAMWARN("You have to set 'HydroBase::timelevels to at least 2");
  }
}
