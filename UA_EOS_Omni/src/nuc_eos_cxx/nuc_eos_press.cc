#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nuc_eos.hh"
#include "helpers.hh"

namespace nuc_eos {

extern "C"
void CCTK_FNAME(nuc_eos_m_kt1_press_eps)(const int *restrict n_in,
			     const double *restrict rho,
			     const double *restrict temp,
			     const double *restrict ye,
			     double *restrict eps,
			     double *restrict prs,
			     int *restrict keyerr,
			     int *restrict anyerr)
{

  using namespace nuc_eos;

  const int n = *n_in;

  *anyerr = 0;
  for(int i=0;i<n;i++) {
    // check if we are fine
    keyerr[i] = checkbounds(rho[i], temp[i], ye[i]);
    if(CCTK_BUILTIN_EXPECT(keyerr[i] != 0, false)) {
      *anyerr = 1;
    }
  }

  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.
  if(*anyerr) return;

  for(int i=0;i<n;i++) {
    int idx[8];
    double delx,dely,delz;
    const double xrho = log(rho[i]);
    const double xtemp = log(temp[i]);
//    const double xrho = rho[i];
//    const double xtemp = temp[i];

    get_interp_spots(xrho,xtemp,ye[i],&delx,&dely,&delz,idx);

    {
      const int iv = 0;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
    }

    {
      const int iv = 1;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(eps[i]),iv);
    }
  }

  // now get rid of ln:
  for(int i=0;i<n;i++) {
    prs[i] = exp(prs[i]);
    eps[i] = exp(eps[i]) - energy_shift;
  }


  return;
}

extern "C"
void CCTK_FNAME(nuc_eos_m_kt0_press)(const int *restrict n_in,
			 const double *restrict rho,
			 double *restrict temp,
			 const double *restrict ye,
			 const double *restrict eps,
			 double *restrict prs,
			 const double *restrict prec,
			 int *restrict keyerr,
			 int *restrict anyerr)
{

  using namespace nuc_eos;

  const int n = *n_in;

  *anyerr = 0;

  for(int i=0;i<n;i++) {
    // check if we are fine
    // Note that this code now requires that the
    // temperature guess be within the table bounds
    keyerr[i] = checkbounds_kt0_noTcheck(rho[i], ye[i]);
    if(CCTK_BUILTIN_EXPECT(keyerr[i] != 0, false)){
      *anyerr = 1;
    }
  }

  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.
  if(*anyerr) return;

  // first must find the temperature
  for(int i=0;i<n;i++) {
    const double lr = log(rho[i]);
    const double lt = log(MIN(MAX(temp[i],eos_tempmin),eos_tempmax));
    double ltout;
    const double epstot = eps[i]+energy_shift;
    if(CCTK_BUILTIN_EXPECT(epstot>0.0e0, true)) {
      // this is the standard scenario; eps is larger than zero
      // and we can operate with logarithmic tables
      const double lxeps = log(epstot);
#if DEBUG
      fprintf(stderr,"%d %15.6E %15.6E %15.6E %15.6E\n",i,lr,lt,ye[i],lxeps);
      fprintf(stderr,"%d %15.6E %15.6E %15.6E %15.6E\n",i,
	      exp(lr),exp(lt),ye[i],exp(lxeps));
#endif
      nuc_eos_findtemp(lr,lt,ye[i],lxeps,*prec,
		       (double *restrict)(&ltout),&keyerr[i]);
    } else {
      // will be overwritten further down, only marks error
      keyerr[i] = 667;
    } // epstot > 0.0

    if(CCTK_BUILTIN_EXPECT(keyerr[i] != 0, false)) {
      // now try negative temperature treatment
      double eps0, eps1;
      int idx[8];
      double delx,dely,delz;

      get_interp_spots_linT_low_eps(lr,temp1,ye[i],&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one_linT_low_eps(idx,delx,dely,delz,&(eps1));

      get_interp_spots_linT_low_eps(lr,temp0,ye[i],&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one_linT_low_eps(idx,delx,dely,delz,&(eps0));

      temp[i] = (epstot-eps0) * (temp1-temp0)/(eps1-eps0) + temp0;

      // set error codes
      *anyerr = 1;
      keyerr[i] = 668;

      // get pressure
      {
	const int iv = 0;
	get_interp_spots_linT_low(lr,temp[i],ye[i],&delx,&dely,&delz,idx);
	nuc_eos_C_linterp_one_linT_low(idx,delx,dely,delz,&(prs[i]),iv);
      }

    } else {
      temp[i] = exp(ltout);
      int idx[8];
      double delx,dely,delz;
      get_interp_spots(lr,ltout,ye[i],&delx,&dely,&delz,idx);
      {
	const int iv = 0;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
      }
    }
  } // loop i<n

  // now get rid of ln:
  for(int i=0;i<n;i++) {
    prs[i] = exp(prs[i]);
  }

  return;
}


extern "C"
void CCTK_FNAME(nuc_eos_P_from_Enthalpy)(const int *restrict n_in,
                                         const double *restrict rho,
                                         double *restrict temp,
                                         const double *restrict ye,
                                         double *restrict enr,
                                         double *restrict enr2,
                                         double *restrict press,
                                         const int *restrict keytemp,
                                         int *restrict keyerr,
                                         const double *restrict prec)
{

  using namespace nuc_eos;

  // keyerr codes:
  // 667 -- no temperature found
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 103 -- temp too high (if keytemp = 1)
  // 104 -- temp too low (if keytemp = 1)
  // 105 -- rho too high
  // 106 -- rho too low

  // keytemp codes:
  // 1 -- coming in with temperature
  // 0 -- coming in with eps, need to find temperature
  // 2 -- coming in with spec. entropy, need to find temperature
  //      (not currently implemented)


  int npoints = 1;
  int anyerr = 0;
  *keyerr = 0;

  // set up local vars

  double eps2 = 0.0;
  double press2 = 0.0;
  double eps = *enr;
  double lr = log(*rho);
  double lt = log(*temp);
  double lenthalpy = eps;

  CCTK_FNAME(nuc_eos_m_kt1_press_eps)(&npoints, rho, temp,ye,&eps2,&press2,keyerr,&anyerr);

  double enthalpy = eps2+press2 / *rho;
  enthalpy = enthalpy + enthalpy*0.01;

  if(*keytemp == 0) {
    double nlt = 0.0;
    nuc_eos_findtemp_enthalpy(lr,lt,*ye,lenthalpy,*prec,&nlt,keyerr);

//    if(*keyerr != 0) {
//      printf("3 Keyerr: %d, Temp: %g Rho: %g Eps: %g Ye: %g\n",*keyerr,*temp, *rho, eps, *ye);
//    }
    lt = nlt;

    *temp = exp(lt);
  } else if(*keytemp == 1) {

  }

  int idx[8];
  double delx,dely,delz,lpress,leps;

  get_interp_spots(lr,lt,*ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&lpress,0);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&leps,1);
  // assign results

  *enr2 = exp(leps) - energy_shift;
  *press= exp(lpress);

  return;

}


} // namespace nuc_eos
