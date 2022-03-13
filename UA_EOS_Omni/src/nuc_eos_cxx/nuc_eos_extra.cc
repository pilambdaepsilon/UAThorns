#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nuc_eos.hh"
#include "helpers.hh"

extern "C" {
  void nuc_eos_get_lims(double *restrict rho_min,
                        double *restrict rho_max,
                        double *restrict t_min,
                        double *restrict t_max,
                        double *restrict ye_min,
                        double *restrict ye_max,
                        double *restrict eps_min,
                        double *restrict eps_max,
                        double *restrict press_min,
                        double *restrict press_max)
  {

    using namespace nuc_eos;


    *rho_min = eos_rhomin;
    *rho_max = eos_rhomax;
    *t_min   = eos_tempmin;
    *t_max   = eos_tempmax;
    *ye_min  = eos_yemin;
    *ye_max  = eos_yemax;
    *eps_min = eos_epsmin;
    *eps_max = eos_epsmax;
    *press_min = eos_pressmin;
    *press_max = eos_pressmax;

    return;
  }
}

extern "C" {

  void CCTK_FNAME (nuc_eos_m_kt1_dpdrhot_dedrhot)(const int *restrict n_in,
                                                  const double *restrict rho,
                                                  double *restrict temp,
                                                  const double *restrict ye,
                                                  double *restrict eps,
                                                  double *restrict eps2,
                                                  double *restrict press,
                                                  double *restrict dpdrho,
                                                  double *restrict dpdt,
                                                  double *restrict dedrho,
                                                  double *restrict dedt,
                                                  const double *restrict prec,
                                                  int *restrict keyerr,
                                                  int *restrict anyerr)
  {

    using namespace nuc_eos;

    const int n = *n_in;
    int keyerr2;

    for(int i=0;i<n;i++) {
      // check if we are fine
      keyerr2 = checkbounds(rho[i], temp[i], ye[i]);
      if(keyerr2 != 0) {
        *anyerr = 1;
        printf("Inside kt1_dpdrhot_dedrhot: keyerr: %d ,[rho]%g, [temp]%g, [ye]%g \n",keyerr2, rho[i], temp[i], ye[i]);
      }
    }

    // Abort if there is any error in checkbounds.
    // This should never happen and the program should abort with
    // a fatal error anyway. No point in doing any further EOS calculations.
    if(*anyerr) return;

    for(int i=0;i<n;i++) {

      int idx[11];
      double delx,dely,delz;
      const double lr = log(rho[i]);
      const double lt = log(MIN(MAX(temp[i],eos_tempmin),eos_tempmax));

      get_interp_spots_d(lr,lt,ye[i],&delx,&dely,&delz,idx);

      // get prs and eps derivatives by rho and temp

      {
        int iv = 0;
        nuc_eos_C_linterp_one_d(idx,delx,dely,delz,iv,&(press[i]), &(dedrho[i]), &(dedt[i]), &(dpdrho[i]), &(dpdt[i]));
        iv = 1;
        nuc_eos_C_linterp_one_d(idx,delx,dely,delz,iv,&(eps[i]), &(dedrho[i]), &(dedt[i]), &(dpdrho[i]), &(dpdt[i]));

        press[i] = exp(press[i]);
        eps2[i] = exp(eps[i]);
        eps[i] = exp(eps[i])-energy_shift;
      }

    }

    return;
  }


}
