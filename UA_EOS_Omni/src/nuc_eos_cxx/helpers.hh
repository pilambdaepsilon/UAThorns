
// helpers.hh edited by Lorenzo Sala

#include "nuc_eos.hh"
#include <cstdlib>
#include <iostream>

namespace nuc_eos {

static inline __attribute__((always_inline))  
int checkbounds(const double xrho, 
		const double xtemp, 
		const double xye) {
  
  using namespace nuc_eos;

  // keyerr codes:
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 103 -- temp too high (if keytemp = 1)
  // 104 -- temp too low (if keytemp = 1)
  // 105 -- rho too high
  // 106 -- rho too low

  if(CCTK_BUILTIN_EXPECT(xrho > eos_rhomax, false)) {
    return 105;
  }
  if(CCTK_BUILTIN_EXPECT(xrho < eos_rhomin, false)) {
    return 106;
  }
  if(CCTK_BUILTIN_EXPECT(xye > eos_yemax, false)) {
    return 101;
  }
  if(CCTK_BUILTIN_EXPECT(xye < eos_yemin, false)) {
    // this is probably not pure and should be removed
    //fprintf(stderr,"xye: %15.6E eos_yemin: %15.6E\n",xye,eos_yemin);
    return 102;
  }
  if(CCTK_BUILTIN_EXPECT(xtemp > eos_tempmax, false)) {
    return 103;
  }
  if(CCTK_BUILTIN_EXPECT(xtemp < eos_tempmin, false)) {
    return 104;
  }
  return 0;
}

static inline __attribute__((always_inline))  
int checkbounds_kt0_noTcheck(const double xrho, 
		const double xye) {
  
  using namespace nuc_eos;

  // keyerr codes:
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 105 -- rho too high
  // 106 -- rho too low

  if(CCTK_BUILTIN_EXPECT(xrho > eos_rhomax, false)) {
    return 105;
  }
  if(CCTK_BUILTIN_EXPECT(xrho < eos_rhomin, false)) {
    return 106;
  }
  if(CCTK_BUILTIN_EXPECT(xye > eos_yemax, false)) {
    return 101;
  }
  if(CCTK_BUILTIN_EXPECT(xye < eos_yemin, false)) {
    return 102;
  }
  return 0;
}


static inline __attribute__((always_inline))  
void get_interp_spots(const double x,
		      const double y,
		      const double z,
		      double* restrict delx,
		      double* restrict dely,
		      double* restrict delz,
		      int* restrict idx) 
{
  using namespace nuc_eos;
  using namespace nuc_eos_private;

  int ix = 1 + (int)( (x - logrho[0] - 1.0e-10) * drhoi );
  int iy = 1 + (int)( (y - logtemp[0] - 1.0e-10) * dtempi );
  int iz = 1 + (int)( (z - yes[0] - 1.0e-10) * dyei );

  ix = MAX( 1, MIN( ix, nrho-1 ) );
  iy = MAX( 1, MIN( iy, ntemp-1 ) );
  iz = MAX( 1, MIN( iz, nye-1 ) );

  idx[0] = NTABLES*(ix + nrho*(iy + ntemp*iz));
  idx[1] = NTABLES*((ix-1) + nrho*(iy + ntemp*iz));
  idx[2] = NTABLES*(ix + nrho*((iy-1) + ntemp*iz));
  idx[3] = NTABLES*(ix + nrho*(iy + ntemp*(iz-1)));
  idx[4] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*iz));
  idx[5] = NTABLES*((ix-1) + nrho*(iy + ntemp*(iz-1)));
  idx[6] = NTABLES*(ix + nrho*((iy-1) + ntemp*(iz-1)));
  idx[7] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*(iz-1)));

   // set up aux vars for interpolation
  *delx = logrho[ix] - x;
  *dely = logtemp[iy] - y;
  *delz = yes[iz] - z;

  return;
}

static inline __attribute__((always_inline))
void get_interp_spots_d(const double x,
		      const double y,
		      const double z,
		      double* restrict delx,
		      double* restrict dely,
		      double* restrict delz,
		      int* restrict idx)
{

  using namespace nuc_eos;
  using namespace nuc_eos_private;

  int ix = 1 + (int)( (x - logrho[0] - 1.0e-10) * drhoi );
  int iy = 1 + (int)( (y - logtemp[0] - 1.0e-10) * dtempi );
  int iz = 1 + (int)( (z - yes[0] - 1.0e-10) * dyei );

  ix = MAX( 1, MIN( ix, nrho-1 ) );
  iy = MAX( 1, MIN( iy, ntemp-1 ) );
  iz = MAX( 1, MIN( iz, nye-1 ) );

  idx[0] = NTABLES*(ix + nrho*(iy + ntemp*iz));
  idx[1] = NTABLES*((ix-1) + nrho*(iy + ntemp*iz));
  idx[2] = NTABLES*(ix + nrho*((iy-1) + ntemp*iz));
  idx[3] = NTABLES*(ix + nrho*(iy + ntemp*(iz-1)));
  idx[4] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*iz));
  idx[5] = NTABLES*((ix-1) + nrho*(iy + ntemp*(iz-1)));
  idx[6] = NTABLES*(ix + nrho*((iy-1) + ntemp*(iz-1)));
  idx[7] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*(iz-1)));
  idx[8] = NTABLES*((ix+1) + nrho*(iy + ntemp*iz));
  idx[9] = NTABLES*(ix + nrho*((iy+1) + ntemp*iz));
  idx[10] = NTABLES*(ix + nrho*(iy + ntemp*(iz+1)));

   // set up aux vars for interpolation
  *delx = logrho[ix] - x;
  *dely = logtemp[iy] - y;
  *delz = yes[iz] - z;

  return;
}

static inline __attribute__((always_inline))  
void get_interp_spots_linT_low(const double x,
			       const double y,
			       const double z,
			       double* restrict delx,
			       double* restrict dely,
			       double* restrict delz,
			       int* restrict idx) 
{
  using namespace nuc_eos;
  using namespace nuc_eos_private;

  int ix = 1 + (int)( (x - logrho[0] - 1.0e-10) * drhoi );
  int iy = 1;
  int iz = 1 + (int)( (z - yes[0] - 1.0e-10) * dyei );

  ix = MAX( 1, MIN( ix, nrho-1 ) );
  iz = MAX( 1, MIN( iz, nye-1 ) );

  idx[0] = NTABLES*(ix + nrho*(iy + ntemp*iz));
  idx[1] = NTABLES*((ix-1) + nrho*(iy + ntemp*iz));
  idx[2] = NTABLES*(ix + nrho*((iy-1) + ntemp*iz));
  idx[3] = NTABLES*(ix + nrho*(iy + ntemp*(iz-1)));
  idx[4] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*iz));
  idx[5] = NTABLES*((ix-1) + nrho*(iy + ntemp*(iz-1)));
  idx[6] = NTABLES*(ix + nrho*((iy-1) + ntemp*(iz-1)));
  idx[7] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*(iz-1)));

   // set up aux vars for interpolation
  *delx = logrho[ix] - x;
  *dely = temp1 - y;
  *delz = yes[iz] - z;

  return;
}

static inline __attribute__((always_inline))
void nuc_eos_C_linterp_one_d(const int* restrict idx,
			   const double delx,
			   const double dely,
			   const double delz,
			   const int iv,
			   double* restrict f,
                           double* restrict depsdrho,
			   double* restrict depsdt,
  			   double* restrict dPdrho,
			   double* restrict dPdt)
{

  using namespace nuc_eos;
  using namespace nuc_eos_private;


  // helper variables
  double fh[11], a[11];

  fh[0] = alltables[iv+idx[0]];
  fh[1] = alltables[iv+idx[1]];
  fh[2] = alltables[iv+idx[2]];
  fh[3] = alltables[iv+idx[3]];
  fh[4] = alltables[iv+idx[4]];
  fh[5] = alltables[iv+idx[5]];
  fh[6] = alltables[iv+idx[6]];
  fh[7] = alltables[iv+idx[7]];
  fh[8] = alltables[iv+idx[8]];
  fh[9] = alltables[iv+idx[9]];
  fh[10] = alltables[iv+idx[10]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = drhoi *   ( fh[1] - fh[0] );
  a[2] = dtempi *   ( fh[2] - fh[0] );
  a[3] = dyei *   ( fh[3] - fh[0] );
  a[4] = drhotempi *  ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = drhoyei *  ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = dtempyei *  ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = drhotempyei * ( fh[7] - fh[0] + fh[1] + fh[2] +
			 fh[3] - fh[4] - fh[5] - fh[6] );
  a[8] = 0.5 * drhoi *   ( fh[1] - fh[8] );
  a[9] = 0.5 * dtempi *   ( fh[2] - fh[9] );

  *f = a[0] + a[1] * delx
    + a[2] * dely
    + a[3] * delz
    + a[4] * delx * dely
    + a[5] * delx * delz
    + a[6] * dely * delz
    + a[7] * delx * dely * delz;

  if (iv==0) {
    *dPdt = -a[2];
    *dPdrho = -a[1];
  }
  if (iv==1) {
    *depsdt = -a[2];
    *depsdrho = -a[1];
  }
  return;
}

static inline __attribute__((always_inline))  
void get_interp_spots_linT_low_eps(const double x,
				   const double y,
				   const double z,
				   double* restrict delx,
				   double* restrict dely,
				   double* restrict delz,
				   int* restrict idx) 
{
  using namespace nuc_eos;
  using namespace nuc_eos_private;

  int ix = 1 + (int)( (x - logrho[0] - 1.0e-10) * drhoi );
  int iy = 1;
  int iz = 1 + (int)( (z - yes[0] - 1.0e-10) * dyei );

  ix = MAX( 1, MIN( ix, nrho-1 ) );
  iz = MAX( 1, MIN( iz, nye-1 ) );

  idx[0] = (ix + nrho*(iy + ntemp*iz));
  idx[1] = ((ix-1) + nrho*(iy + ntemp*iz));
  idx[2] = (ix + nrho*((iy-1) + ntemp*iz));
  idx[3] = (ix + nrho*(iy + ntemp*(iz-1)));
  idx[4] = ((ix-1) + nrho*((iy-1) + ntemp*iz));
  idx[5] = ((ix-1) + nrho*(iy + ntemp*(iz-1)));
  idx[6] = (ix + nrho*((iy-1) + ntemp*(iz-1)));
  idx[7] = ((ix-1) + nrho*((iy-1) + ntemp*(iz-1)));

   // set up aux vars for interpolation
  *delx = logrho[ix] - x;
  *dely = temp1 - y;
  *delz = yes[iz] - z;

  return;
}


static inline __attribute__((always_inline))  
void nuc_eos_C_linterp_one(const int* restrict idx, 
			   const double delx, 
			   const double dely,
			   const double delz,
			   double* restrict f,
			   const int iv)
{
  using namespace nuc_eos;
  using namespace nuc_eos_private;

  // helper variables
  double fh[8], a[8];

  fh[0] = alltables[iv+idx[0]];
  fh[1] = alltables[iv+idx[1]];
  fh[2] = alltables[iv+idx[2]];
  fh[3] = alltables[iv+idx[3]];
  fh[4] = alltables[iv+idx[4]];
  fh[5] = alltables[iv+idx[5]];
  fh[6] = alltables[iv+idx[6]];
  fh[7] = alltables[iv+idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = drhoi *   ( fh[1] - fh[0] );
  a[2] = dtempi *   ( fh[2] - fh[0] );
  a[3] = dyei *   ( fh[3] - fh[0] );
  a[4] = drhotempi *  ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = drhoyei *  ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = dtempyei *  ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = drhotempyei * ( fh[7] - fh[0] + fh[1] + fh[2] + 
			 fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0] + a[1] * delx
    + a[2] * dely
    + a[3] * delz
    + a[4] * delx * dely
    + a[5] * delx * delz
    + a[6] * dely * delz
    + a[7] * delx * dely * delz;

  return;
}

static inline __attribute__((always_inline))  
void nuc_eos_C_linterp_one_linT_low(const int* restrict idx, 
				    const double delx, 
				    const double dely,
				    const double delz,
				    double* restrict f,
				    const int iv)
{
  using namespace nuc_eos;
  using namespace nuc_eos_private;

  // helper variables
  double fh[8], a[8];

  fh[0] = alltables[iv+idx[0]];
  fh[1] = alltables[iv+idx[1]];
  fh[2] = alltables[iv+idx[2]];
  fh[3] = alltables[iv+idx[3]];
  fh[4] = alltables[iv+idx[4]];
  fh[5] = alltables[iv+idx[5]];
  fh[6] = alltables[iv+idx[6]];
  fh[7] = alltables[iv+idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = drhoi *   ( fh[1] - fh[0] );
  a[2] = dlintempi *   ( fh[2] - fh[0] );
  a[3] = dyei *   ( fh[3] - fh[0] );
  a[4] = drholintempi *  ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = drhoyei *  ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = dlintempyei *  ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = drholintempyei * ( fh[7] - fh[0] + fh[1] + fh[2] + 
			 fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0] + a[1] * delx
    + a[2] * dely
    + a[3] * delz
    + a[4] * delx * dely
    + a[5] * delx * delz
    + a[6] * dely * delz
    + a[7] * delx * dely * delz;

  return;
}


static inline __attribute__((always_inline))  
void nuc_eos_C_linterp_one_linT_low_eps(const int* restrict idx, 
					const double delx, 
					const double dely,
					const double delz,
					double* restrict f)

{
  using namespace nuc_eos;
  using namespace nuc_eos_private;

  // helper variables
  double fh[8], a[8];

  fh[0] = epstable[idx[0]];
  fh[1] = epstable[idx[1]];
  fh[2] = epstable[idx[2]];
  fh[3] = epstable[idx[3]];
  fh[4] = epstable[idx[4]];
  fh[5] = epstable[idx[5]];
  fh[6] = epstable[idx[6]];
  fh[7] = epstable[idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = drhoi *   ( fh[1] - fh[0] );
  a[2] = dlintempi *   ( fh[2] - fh[0] );
  a[3] = dyei *   ( fh[3] - fh[0] );
  a[4] = drholintempi *  ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = drhoyei *  ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = dlintempyei *  ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = drholintempyei * ( fh[7] - fh[0] + fh[1] + fh[2] + 
			 fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0] + a[1] * delx
    + a[2] * dely
    + a[3] * delz
    + a[4] * delx * dely
    + a[5] * delx * delz
    + a[6] * dely * delz
    + a[7] * delx * dely * delz;

  return;
}




static inline __attribute__((always_inline))
double linterp2D(const double *restrict xs, 
		 const double *restrict ys, 
		 const double *restrict fs, 
		 const double x, 
		 const double y)
{

  //  2     3 
  //
  //  0     1
  //
  // first interpolate in x between 0 and 1, 2 and 3
  // then interpolate in y
  // assume rectangular grid
  
  double dxi = 1./(xs[1]-xs[0]);
  double dyi = 1./(ys[1]-ys[0]); // x*1./y uses faster instructions than x/y
  double t1 = (fs[1]-fs[0])*dxi * (x - xs[0]) + fs[0];
  double t2 = (fs[3]-fs[2])*dxi * (x - xs[0]) + fs[2];

  return (t2 - t1)*dyi * (y-ys[0]) + t1;
}

static inline __attribute__((always_inline))
void bisection(const double lr, 
	       const double lt0,
	       const double ye,
	       const double leps0,
	       const double prec,
	       double *restrict ltout,
	       const int iv,
	       int *restrict keyerrt) {
  // iv is the index of the variable we do the bisection on

  using namespace nuc_eos;
  using namespace nuc_eos_private;

  int bcount = 0; 
  int maxbcount = 80;
  int itmax = 50;

  const double dlt0p = log(1.1);
  const double dlt0m = log(0.9);
  const double dltp = log(1.2);
  const double dltm = log(0.8);

  double leps0_prec = fabs(leps0*prec);

  // temporary local vars
  double lt, lt1, lt2;
  double ltmin = logtemp[0];
  double ltmax = logtemp[ntemp-1];
  double f1,f2,fmid,dlt,ltmid;
  double f1a = 0.0;
  double f2a = 0.0;
  double delx,dely,delz;
  int idx[8];
  
  // LSMOD (Modification made by Lorenzo Sala)
  // LSMOD: The following lines calculate eps in 
  //        f2a = eps(rho,Tmin, Ye) and f1a = eps(rho,Tmax,Ye)
  get_interp_spots(lr,ltmax,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a,iv);
  get_interp_spots(lr,ltmin,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a,iv);

  // prepare
  // check if your energy is actually tabulated at this rho and ye.
  // f2a is the energy evaluated at ltmin, so it is the minimum energy tabulated 
  // at this rho ad ye.
  // If leps0 <= f2a, then ltout is likely to be the minimum temperature tabulated.
  if(leps0 <= f2a) { // + 1.0E-6
    *ltout = ltmin;
    return;
  }

  /* // If leps0 >= f1a, then ltout is likely to be the maximum temperature tabulated.
  if(leps0 >= f1a) { // + 1.0E-6
    *ltout = ltmax;
    return;
  } */

  // otherwise, proceed finding extrema for applying bisection method.
  lt = lt0;
  lt1 = MIN(lt0 + dlt0p,ltmax);
  lt2 = MAX(lt0 + dlt0m,ltmin);

  get_interp_spots(lr,lt1,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a,iv);  

  get_interp_spots(lr,lt2,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a,iv);  

  f1=f1a-leps0;
  f2=f2a-leps0;

  // iterate until we bracket the right eps, but enforce
  // dE/dt > 0, so eps(lt1) > eps(lt2)
  while(f1*f2 >= 0.0) {
    lt1 = MIN(lt1 + dltp,ltmax);
    lt2 = MAX(lt2 + dltm,ltmin);
    get_interp_spots(lr,lt1,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a,iv);  
    
    get_interp_spots(lr,lt2,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a,iv);  

    f1=f1a-leps0;
    f2=f2a-leps0;

#if DEBUG
    fprintf(stderr,"bisection bracketing it %d, f1: %15.6E, f2: %15.6E, lt1: %15.6E, lt2: %15.6E, f1a: %18.11E, f2a: %18.11E leps0: %18.11E\n",
	    bcount,f1,f2,lt1,lt2,f1a,f2a,leps0);
#endif

    bcount++;
    if(CCTK_BUILTIN_EXPECT(bcount >= maxbcount, false)) {
#if DEBUG
      fprintf(stderr,"bcount out of range it %d, lr: %15.6E, lt1: %15.6E, lt2: %15.6E, f1a: %18.11E, f2a: %18.11E leps0: %18.11E, ye: %15.6E\n",
	      bcount,lr,lt1,lt2,f1a,f2a,leps0,ye);
#endif
      *keyerrt = 667;
      return;
    }
  } // while

  if(f1 < 0.0) {
    lt = lt1;
    dlt = lt2 - lt1;
  } else {
    lt = lt2;
    dlt = lt1 - lt2;
  }

#if DEBUG
    fprintf(stderr,"bisection step 2 it -1, fmid: %15.6E ltmid: %15.6E dlt: %15.6E\n",
	    f2,lt,dlt);
    fprintf(stderr,"ltmax: %15.6E\n",ltmax);
#endif

  int it;
  for(it=0;it<itmax;it++) {
    dlt = dlt * 0.5;
    ltmid = lt + dlt;
    get_interp_spots(lr,ltmid,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a,iv); 
    
    fmid=f2a-leps0;
    if(CCTK_BUILTIN_EXPECT(fmid <= 0.0, false)) lt=ltmid;
#if DEBUG
    fprintf(stderr,"bisection step 2 it %d, fmid: %15.6E f2a: %15.6E lt: %15.6E ltmid: %15.6E dlt: %15.6E\n",
	    it,fmid,f2a,lt,ltmid,dlt);
#endif

    if(CCTK_BUILTIN_EXPECT(fabs(leps0-f2a) <= leps0_prec, false)) {
      *ltout = ltmid;
      return;
    }
  } // for it = 0

  *keyerrt = 667;
  return;
} // bisection



static inline __attribute__((always_inline))
void nuc_eos_findtemp(const double lr, 
		      const double lt0,
		      const double ye,
		      const double lepsin,
		      const double prec,
		      double *restrict ltout,
		      int *keyerrt) {

  using namespace nuc_eos;
  using namespace nuc_eos_private;

  // local variables
  const int itmax = 200; // use at most 10 iterations, then go to bisection
  double dlepsdlti; // 1 / derivative dlogeps/dlogT
  double ldt;
  double leps,leps0; // temp vars for eps
  double ltn, lt; // temp vars for temperature
  double oerr;
  const double ltmax = logtemp[ntemp-1]; // max temp
  const double ltmin = logtemp[0]; // min temp
  const int iv = 1;
  int it = 0;

  // setting up some vars
  *keyerrt = 0;
  leps0 = lepsin;
  lt = lt0;

  // step 1: do we already have the right temperature
  int idx[8];
  double delx,dely,delz;
  get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&leps,iv);
#if DEBUG
  fprintf(stderr,"it: %d t: %15.6E leps: %15.6E eps0: %15.6E del: %15.6E\n",
	  it,lt, leps,leps0,fabs(leps-leps0)/(fabs(leps0)));
#endif
  // TODO: profile this to see which outcome is more likely
  if(fabs(leps-leps0) < prec*fabs(leps0)) {
    *ltout = lt0;
    return;
  }

  oerr = 1.0e90;
  double fac = 1.0;
  const int irho = MIN(MAX(1 + (int)(( lr - logrho[0] - 1.0e-12) * drhoi),1),nrho-1); 
  const int iye = MIN(MAX(1 + (int)(( ye - yes[0] - 1.0e-12) * dyei),1),nye-1); 

  /* ******* if temp low for high density, switch directly to bisection. 
             Verifying Newton-Raphson result evaluating the derivative. 
             The variable shouldgotobisection will be modified accordingly
             to the value of derivative of eps wrt temp ******* */
  bool shouldgotobisection = false; // LSMOD
  while(it < itmax && shouldgotobisection == false) {
    it++;

    // step 2: check if the two bounding values of the temperature
    //         give eps values that enclose the new eps.
    const int itemp = MIN(MAX(1 + (int)(( lt - logtemp[0] - 1.0e-12) * dtempi),1),ntemp-1); 

    double epst1, epst2;
    // lower temperature
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = 1 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye-1)));
      fs[0] = alltables[ifs];
      // point 1 
      ifs = 1 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye-1)));
      fs[1] = alltables[ifs];
      // point 2 
      ifs = 1 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye)));
      fs[2] = alltables[ifs];
      // point 3
      ifs = 1 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye)));
      fs[3] = alltables[ifs];
      
      epst1 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
    }
    // upper temperature
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = 1 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye-1)));
      fs[0] = alltables[ifs];
      // point 1 
      ifs = 1 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye-1)));
      fs[1] = alltables[ifs];
      // point 2 
      ifs = 1 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye)));
      fs[2] = alltables[ifs];
      // point 3
      ifs = 1 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye)));
      fs[3] = alltables[ifs];
      
      epst2 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
    }

    // Check if we are already bracketing the input internal
    // energy. If so, interpolate for new T.
    if(CCTK_BUILTIN_EXPECT((leps0 - epst1) * (leps0 - epst2) <= 0., false)) {
      
      *ltout = (logtemp[itemp]-logtemp[itemp-1]) / (epst2 - epst1) * 
	(leps0 - epst1) + logtemp[itemp-1];
     
      return;
    }

    // well, then do a Newton-Raphson step
    // first, guess the derivative
    dlepsdlti = (logtemp[itemp]-logtemp[itemp-1])/(epst2-epst1);
    ldt = -(leps - leps0) * dlepsdlti * fac;

    //LSMOD: too large a dlt means that the energy dependence on the temperature
    //       is weak ==> We'd better try bisection.
    //       Factor 1/12.0 come from tests by LSMOD
    //       This is done in order to limit the "velocity" of T variation
    //       given by Newton-Raphson.    
    if(ldt > (ltmax-ltmin) / 12.0 ) shouldgotobisection = true;

    ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
    lt = ltn;

    get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&leps,iv);

#if DEBUG
    fprintf(stderr,"%d %d %d\n",irho,itemp,iye);
    fprintf(stderr,"it: %d t: %15.6E leps: %15.6E eps0: %15.6E del: %15.6E\n",
	    it,lt, leps,leps0,fabs(leps-leps0)/(fabs(leps0)));
#endif
    // drive the thing into the right direction
    double err = fabs(leps-leps0);
    if(oerr < err) fac *= 0.9;
    oerr = err;

    if(CCTK_BUILTIN_EXPECT(err < prec*fabs(leps0), false)) {
      *ltout = lt;
      return;
    }



  } // while(it < itmax)

  // try bisection
  // bisection(lr, lt0, ye, leps0, ltout, 1, prec, keyerrt);
#if DEBUG
  fprintf(stderr, "Failed to converge. This is bad. Trying bisection!\n");
#endif
  bisection(lr,lt0,ye,leps0,prec,ltout,1,keyerrt);
#if DEBUG
  if(*keyerrt==667) {
     fprintf(stderr,"This is worse. Bisection failed!\n");
     abort();
  }
#endif


  return;
}

static inline __attribute__((always_inline))
void bisection_enthalpy(const double lr,
	       const double lt0,
	       const double ye,
	       double leps0,
	       const double prec,
	       double *restrict ltout,
	       const int iv,
	       int *restrict keyerrt) {
  // iv is the index of the variable we do the bisection on

  using namespace nuc_eos;
  using namespace nuc_eos_private;


  int bcount = 0;
  int bcount2 = 0;
  int maxbcount = 160;
  int itmax = 200;

  const double dlt0p = log(1.1);
  const double dlt0m = log(0.9);
  const double dltp = log(1.2);
  const double dltm = log(0.8);

  const double leps0_prec = leps0*prec;

  // temporary local vars
  double lt, lt1, lt2;
  double ltmin = logtemp[0];
//  double ltmin = log(0.1);
//  double ltmin = log(0.5);
  double ltmax = logtemp[ntemp-1];
  double f1,f2,fmid,dlt,ltmid,fmida_ent;
  double f1a[2] = {0};
  double f2a[2] = {0};
  double fmida[2] = {0};
  double delx,dely,delz;
  int idx[8];

//  energy_shift = 0.0;
  
  // LSMOD (Modification made by Lorenzo Sala)
  // LSMOD: The following lines calculate eps in 
  //        f2a = eps(rho,Tmin, Ye) and f1a = eps(rho,Tmax,Ye)
  get_interp_spots(lr,ltmax,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[0],0);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[1],1);

  get_interp_spots(lr,ltmin,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a[0],0);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a[1],1);

  // prepare
  // check if your energy is actually tabulated at this rho and ye.
  // f2a[0]([1]) is the pressure (energy) evaluated at ltmin, so it is the minimum energy tabulated 
  // at this rho ad ye.
  // If leps0 <= f2a, then ltout is likely to be the minimum temperature tabulated.

  double f1a_ent0 = exp(f1a[1])-energy_shift+exp(f1a[0])/exp(lr);// + 1.0;
  double f2a_ent0 = exp(f2a[1])-energy_shift+exp(f2a[0])/exp(lr);// + 1.0;
  if(leps0 <= f2a_ent0 + 1.0e-6) { // + 1.0E-6
    *ltout = ltmin;
    return;
  }

  /* // If leps0 >= f1a, then ltout is likely to be the maximum temperature tabulated.
  if(leps0 >= f1a) { // + 1.0E-6
    *ltout = ltmax;
    return;
  } */

  // otherwise, proceed finding extrema for applying bisection method.

  // prepare
  lt = lt0;
  lt1 = MIN(lt0 + dlt0p,ltmax);
  lt2 = MAX(lt0 + dlt0m,ltmin);

  get_interp_spots(lr,lt1,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[0],0);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[1],1);

  get_interp_spots(lr,lt2,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a[0],0);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a[1],1);

  double f1a_ent = exp(f1a[1])-energy_shift+exp(f1a[0])/exp(lr);
  double f2a_ent = exp(f2a[1])-energy_shift+exp(f2a[0])/exp(lr);
  f1=f1a_ent-leps0;
  f2=f2a_ent-leps0;
#if DEBUG
  printf("Bisection start: Count = %d, press = %g, eps = %g, f1a_ent = %g, f2a_ent = %g, enthalpy_in = %g\n", bcount,f1a[0], f1a[1], f1a_ent, f2a_ent, leps0);
  printf("Count = %d, f1a_ent = %g, f2a_ent = %g, f1 = %g, f2 = %g, ltmin=%g\n", bcount,f1a_ent,f2a_ent,f1,f2,ltmin);
#endif
  // iterate until we bracket the right eps, but enforce
  // dE/dt > 0, so eps(lt1) > eps(lt2)
  while(f1*f2 >= 0.0) {
    lt1 = MIN(lt1 + dltp,ltmax);
    lt2 = MAX(lt2 + dltm,ltmin);

    get_interp_spots(lr,lt1,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[0],0);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[1],1);

    get_interp_spots(lr,lt2,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a[0],0);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a[1],1);

    f1a_ent = exp(f1a[1])-energy_shift+exp(f1a[0])/exp(lr);
    f2a_ent = exp(f2a[1])-energy_shift+exp(f2a[0])/exp(lr);
    f1=f1a_ent-leps0;
    f2=f2a_ent-leps0;

#if DEBUG
    printf("Count = %d, press = %g, eps = %g, enthalpy_in = %g\n", bcount,f1a[0], f1a[1], leps0);
    printf("Count = %d, press2 = %g, eps2 = %g, enthalpy_in = %g\n", bcount,f2a[0], f2a[1], leps0);
    printf("Count = %d, lt1 = %g, lt2 = %g\n", bcount,lt1,lt2);
    printf("Count = %d, f1a_ent = %g, f2a_ent = %g, f1 = %g, f2 = %g\n", bcount,f1a_ent,f2a_ent,f1,f2);
#endif

#if DEBUG
    fprintf(stderr,"bisection bracketing it %d, f1: %15.6E, f2: %15.6E, lt1: %15.6E, lt2: %15.6E, f1a: %18.11E, f2a: %18.11E leps0: %18.11E\n",
	    bcount,f1,f2,lt1,lt2,f1a_ent,f2a_ent,leps0);
#endif
    bcount++;
    if(bcount >= maxbcount && bcount2 == 0) {
      bcount2 = 1;
      bcount = 0;
      get_interp_spots(lr,ltmin,ye,&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[0],0);
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[1],1);

      double ent_min = exp(f1a[1])-energy_shift+exp(f1a[0])/exp(lr);
      leps0 = fmax(ent_min*0.99,ent_min*1.01);
#if DEBUG
      printf("Adjusting ent_min: Old ent = %g\n", leps0);
      printf("Adjusting ent_min: New ent = %g\n", leps0);
#endif
    }
    else if(bcount >= maxbcount && bcount2 == 1) {
      *keyerrt = 667;
      return;
    }
    else {
    }

  } // while

#if DEBUG
    fprintf(stderr,"bisection step 2 it -1, fmid: %15.6E ltmid: %15.6E dlt: %15.6E\n",
	    f2,lt,dlt);
    fprintf(stderr,"ltmax: %15.6E\n",ltmax);
#endif

  get_interp_spots(lr,lt1,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[0],0);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[1],1);
  f1a_ent = exp(f1a[1])-energy_shift+exp(f1a[0])/exp(lr);

  get_interp_spots(lr,lt2,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a[0],0);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a[1],1);

  f2a_ent = exp(f2a[1])-energy_shift+exp(f2a[0])/exp(lr);
  int it;
  for(it=0;it<itmax;it++) {
    ltmid = log((exp(lt1) + exp(lt2))*0.5);

    get_interp_spots(lr,ltmid,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&fmida[0],0);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&fmida[1],1);
    fmida_ent = exp(fmida[1])-energy_shift+exp(fmida[0])/exp(lr);

    if ((f1a_ent - leps0)*(fmida_ent-leps0) < 0.0) {
	    f2a_ent = fmida_ent;
	    lt2 = ltmid;
    }
    else if ((f2a_ent - leps0)*(fmida_ent-leps0) < 0.0) {
	    f1a_ent = fmida_ent;
	    lt1 = ltmid;
    }
    else {
#if DEBUG
      printf("Bracketing failed for bisection interval\n");
#endif
    }

#if DEBUG
    fprintf(stderr,"bisection step 2 it %d, fmid: %15.6E f2a: %15.6E lt: %15.6E ltmid: %15.6E dlt: %15.6E, crit = %15.6E,leps0_prec = %15.6E\n",
	    it,fmid,f2a_ent,lt,ltmid,dlt,fabs(leps0-f2a_ent),leps0_prec);
#endif

#if DEBUG
     printf("It = %d, press2 = %g, eps2 = %g, enthalpy = %g, enthalpy_in = %g\n", it,f2a[0], f2a[1], f2a_ent, leps0);
     printf("It = %d, ltmid = %g, lt0 = %g\n", it,ltmid,lt0);
     printf("It = %d, f2a_ent = %g, err = %g\n", it,f2a_ent,fabs(1.0-f2a_ent/leps0));
#endif

    if(fabs(leps0-fmida_ent) <= fabs(leps0_prec*5.0e-2)) {

#if DEBUG
      printf("Bisection found correct temp: It = %d, f2a_ent = %g, err = %g, threshold = %g\n", it,f2a_ent,fabs(leps0-f2a_ent),fabs(leps0_prec*1.0e-2));
#endif
      *ltout = ltmid;
      return;
    }
  } // for it = 0
  *keyerrt = 667;
  return;
} 

// bisection
  static inline __attribute__((always_inline))
void nuc_eos_findtemp_enthalpy(const double lr,
		      const double lt0,
		      const double ye,
		      const double enthalpyin,
		      const double prec,
		      double *restrict ltout,
		      int *keyerrt) {

  using namespace nuc_eos;
  using namespace nuc_eos_private;

    // local variables
    const int itmax = 20; // use at most 10 iterations, then go to bisection
    double dentdlti; // 1 / derivative dentropy/dlogT
    double ldt;
    double ent,ent0; // temp vars for eps
    double ltn, lt, lt1; // temp vars for temperature
    double oerr;
    const double ltmax = logtemp[ntemp-1]; // max temp
    const double ltmin = logtemp[0]; // min temp
//    const double ltmin = log(0.1); // min temp
    const int iv = 1;
    int it = 0;

    double dlenthalpydlt,denthalpydtc,dlenthalpydltc,dlpressdlt,dlepsdlt,depsdrho,dpdrho; // derivative dlogenthalpy/dlogT
    double lenthalpy,enthalpy,lenthalpyc,enthalpyc,lenthalpy0,enthalpy0,lenthalpy1,enthalpy1;
    double leps, lpress;
    double presst1, presst2, epst1, epst2, enthalpyt1,enthalpyt2,enthalpyc1,enthalpyc2;
    // setting up some vars
    *keyerrt = 0;
    lenthalpy0 = enthalpyin;
    lenthalpy1 = lenthalpy0;
    lt = lt0;
    lt1 = lt;

    // step 1: do we already have the right temperature
    int idx[11];
    double delx,dely,delz;

    get_interp_spots_d(lr,lt,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one_d(idx,delx,dely,delz,1,&leps,&depsdrho,&dlepsdlt,&dpdrho,&dlepsdlt);
    nuc_eos_C_linterp_one_d(idx,delx,dely,delz,0,&lpress,&depsdrho,&dlepsdlt,&dpdrho,&dlpressdlt);

#if DEBUG
    fprintf(stderr,"it: %d t: %15.6E leps: %15.6E eps0: %15.6E del: %15.6E\n",
	    it,lt,ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
#endif

    enthalpyc = exp(leps) - energy_shift + exp(lpress)/exp(lr);
    lenthalpy = enthalpyc;
    denthalpydtc = exp(leps)/exp(lt)*dlepsdlt + exp(lpress)/exp(lt)*dlpressdlt/exp(lr);
    dlenthalpydlt = exp(lt)/enthalpyc*denthalpydtc;


#if DEBUG
    printf("! eps = %16.15e\n",exp(leps)-energy_shift);
    printf("! press = %16.15e\n",exp(lpress));
    printf("! leps = %16.15e\n",leps);
    printf("! lpress = %16.15e\n",lpress);
    printf("! lr = %16.15e\n",lr);
    printf("! lt = %16.15e\n",lt);
    printf("! r = %16.15e\n",exp(lr));
    printf("! t = %16.15e\n",exp(lt));
    printf("! lenthalpy0 = %16.15e\n",lenthalpy0);
    printf("! lenthalpy1 = %16.15e\n",lenthalpy1);
    printf("! lenthalpy = %16.15e\n",lenthalpy);
    printf("! enthalpyc = %16.15e\n",enthalpyc);
    printf("! dlepsdlt = %16.15e\n",dlepsdlt);
    printf("! dlpressdlt = %16.15e\n",dlpressdlt);
    printf("! dlenthalpydlt = %16.15e\n",dlenthalpydlt);
    printf("! denthalpydtc = %16.15e\n",denthalpydtc);
    printf("------------------------------------\n");
#endif

    if(fabs(lenthalpy-lenthalpy0) < prec*1.0e-2*fabs(lenthalpy0)) {
      *ltout = lt0;
#if DEBUG
      printf("Already have correct temp\n");
#endif
      return;
    }

    lt1 = lt;
    lenthalpy1 = lenthalpy;

    oerr = 1.0e90;
    double fac = 1.0;
    const int irho = MIN(MAX(1 + (int)(( lr - logrho[0] - 1.0e-12) * drhoi),1),nrho-1);
    const int iye = MIN(MAX(1 + (int)(( ye - yes[0] - 1.0e-12) * dyei),1),nye-1);

    while(it < itmax) {
      it++;

      // step 2: check if the two bounding values of the temperature
      //         give eps values that enclose the new eps.
      const int itemp = MIN(MAX(1 + (int)(( lt - logtemp[0] - 1.0e-12) * dtempi),1),ntemp-1);

      double enthalpy1, enthalpy2;
      // lower temperature
      {
	// get data at 4 points
	double fs[4],fs1[4];
	// point 0
	int ifs = 0 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye-1)));
	int ifs1 = 1 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye-1)));
	fs[0] = alltables[ifs];
	fs1[0] = alltables[ifs1];
	// point 1
	ifs = 0 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye-1)));
	ifs1 = 1 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye-1)));
	fs[1] = alltables[ifs];
	fs1[1] = alltables[ifs1];
	// point 2
	ifs = 0 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye)));
	ifs1 = 1 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye)));
	fs[2] = alltables[ifs];
	fs1[2] = alltables[ifs1];
	// point 3
	ifs = 0 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye)));
	ifs1 = 1 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye)));
	fs[3] = alltables[ifs];
	fs1[3] = alltables[ifs1];

	presst1 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
 	epst1 = linterp2D(&logrho[irho-1],&yes[iye-1], fs1, lr, ye);
        enthalpyc1 =  exp(epst1) - energy_shift + exp(presst1)/exp(lr);


      }
      // upper temperature
      {
	// get data at 4 points
	double fs[4],fs1[4];
	// point 0
	int ifs = 0 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye-1)));
	int ifs1 = 1 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye-1)));
	fs[0] = alltables[ifs];
	fs1[0] = alltables[ifs1];
	// point 1
	ifs = 0 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye-1)));
	ifs1 = 1 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye-1)));
	fs[1] = alltables[ifs];
	fs1[1] = alltables[ifs1];
	// point 2
	ifs = 0 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye)));
	ifs = 1 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye)));
	fs[2] = alltables[ifs];
	fs1[2] = alltables[ifs1];
	// point 3
	ifs = 0 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye)));
	ifs1 = 1 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye)));
	fs[3] = alltables[ifs];
	fs1[3] = alltables[ifs1];

 	presst2 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
 	epst2 = linterp2D(&logrho[irho-1],&yes[iye-1], fs1, lr, ye);
        enthalpyc2 = exp(epst2) - energy_shift + exp(presst2)/exp(lr);
      }

      // Check if we are already bracketing the input internal
      // energy. If so, interpolate for new T.
//      if((lenthalpy0 - enthalpyc1) * (lenthalpy0 - enthalpyc2) <= 0.) {
//
//	*ltout = (logtemp[itemp]-logtemp[itemp-1]) / (enthalpyc2 - enthalpyc1) *
//	  (lenthalpy0 - enthalpyc1) + logtemp[itemp-1];
////#if DEBUG
////        printf("Found correct temp in bracketing\n", it);
////#endif
//	return;
//      }

      // well, then do a Newton-Raphson step
      // first, guess the derivative
      denthalpydtc = exp(leps)/exp(lt)*dlepsdlt + exp(lpress)/exp(lt)*dlpressdlt/exp(lr);
      dlenthalpydlt = exp(lt)*denthalpydtc;

      ldt = -(lenthalpy - lenthalpy0) / dlenthalpydlt;

      ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
      lt = ltn;

      get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&leps,1);
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&lpress,0);
      enthalpyc = exp(leps) - energy_shift + exp(lpress)/exp(lr);
      lenthalpy = enthalpyc;


#if DEBUG
      fprintf(stderr,"%d %d %d\n",irho,itemp,iye);
      fprintf(stderr,"it: %d t: %15.6E ent: %15.6E ent0: %15.6E del: %15.6E\n",
	      it,lt, ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
#endif
      // drive the thing into the right direction
      double err = fabs(lenthalpy-lenthalpy0);
      if(oerr < err) fac *= 0.9;
      oerr = err;

      if(err < prec*1.0e-2*fabs(lenthalpy0)) {
	*ltout = lt;
#if DEBUG
        printf("Found correct temp in %d iterations\n", it);
#endif
	return;
      }
    } // while(it < itmax)

    // try bisection
#if DEBUG
    fprintf(stderr, "Failed to converge. This is bad. Trying bisection!\n");
#endif
    bisection_enthalpy(lr,lt0,ye,lenthalpy0,prec,ltout,1,keyerrt);
#if DEBUG
    if(*keyerrt==667) {
      fprintf(stderr,"This is worse. Bisection failed!\n");
    }
#endif


    return;
}

  static inline __attribute__((always_inline))
    void nuc_eos_findtemp_entropy(const double lr, 
				  const double lt0,
				  const double ye,
				  const double entin,
				  const double prec,
				  double *restrict ltout,
				  int *keyerrt) {

    using namespace nuc_eos;
    using namespace nuc_eos_private;

    // local variables
    const int itmax = 200; // use at most 10 iterations, then go to bisection
    double dentdlti; // 1 / derivative dentropy/dlogT
    double ldt;
    double ent,ent0; // temp vars for eps
    double ltn, lt; // temp vars for temperature
    double oerr;
    const double ltmax = logtemp[ntemp-1]; // max temp
    const double ltmin = logtemp[0]; // min temp
    const int iv = 2;
    int it = 0;

    // setting up some vars
    *keyerrt = 0;
    ent0 = entin;
    lt = lt0;

    // step 1: do we already have the right temperature
    int idx[8];
    double delx,dely,delz;
    get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&ent,iv);
#if DEBUG
    fprintf(stderr,"it: %d t: %15.6E leps: %15.6E eps0: %15.6E del: %15.6E\n",
	    it,lt,ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
#endif
    if(CCTK_BUILTIN_EXPECT(fabs(ent-ent0) < prec*fabs(ent0), false)) {
      *ltout = lt0;
      return;
    }

    oerr = 1.0e90;
    double fac = 1.0;
    const int irho = MIN(MAX(1 + (int)(( lr - logrho[0] - 1.0e-12) * drhoi),1),nrho-1); 
    const int iye = MIN(MAX(1 + (int)(( ye - yes[0] - 1.0e-12) * dyei),1),nye-1); 

    while(it < itmax) {
      it++;

      // step 2: check if the two bounding values of the temperature
      //         give eps values that enclose the new eps.
      const int itemp = MIN(MAX(1 + (int)(( lt - logtemp[0] - 1.0e-12) * dtempi),1),ntemp-1); 

      double ent1, ent2;
      // lower temperature
      {
	// get data at 4 points
	double fs[4];
	// point 0
	int ifs = 2 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye-1)));
	fs[0] = alltables[ifs];
	// point 1 
	ifs = 2 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye-1)));
	fs[1] = alltables[ifs];
	// point 2 
	ifs = 2 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye)));
	fs[2] = alltables[ifs];
	// point 3
	ifs = 2 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye)));
	fs[3] = alltables[ifs];
      
	ent1 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
      }
      // upper temperature
      {
	// get data at 4 points
	double fs[4];
	// point 0
	int ifs = 2 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye-1)));
	fs[0] = alltables[ifs];
	// point 1 
	ifs = 2 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye-1)));
	fs[1] = alltables[ifs];
	// point 2 
	ifs = 2 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye)));
	fs[2] = alltables[ifs];
	// point 3
	ifs = 2 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye)));
	fs[3] = alltables[ifs];
      
	ent2 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
      }

      // Check if we are already bracketing the input internal
      // energy. If so, interpolate for new T.
      if(CCTK_BUILTIN_EXPECT((ent0 - ent1) * (ent0 - ent2) <= 0., false)) {
      
	*ltout = (logtemp[itemp]-logtemp[itemp-1]) / (ent2 - ent1) * 
	  (ent0 - ent1) + logtemp[itemp-1];
     
	return;
      }

      // well, then do a Newton-Raphson step
      // first, guess the derivative
      dentdlti = (logtemp[itemp]-logtemp[itemp-1])/(ent2-ent1);
      ldt = -(ent - ent0) * dentdlti * fac;

      ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
      lt = ltn;

      get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&ent,iv);

#if DEBUG
      fprintf(stderr,"%d %d %d\n",irho,itemp,iye);
      fprintf(stderr,"it: %d t: %15.6E ent: %15.6E ent0: %15.6E del: %15.6E\n",
	      it,lt, ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
#endif
      // drive the thing into the right direction
      double err = fabs(ent-ent0);
      if(oerr < err) fac *= 0.9;
      oerr = err;

      if(CCTK_BUILTIN_EXPECT(err < prec*fabs(ent0), false)) {
	*ltout = lt;
	return;
      }



    } // while(it < itmax)

    // try bisection
#if DEBUG
    fprintf(stderr, "Failed to converge. This is bad. Trying bisection!\n");
#endif
    bisection(lr,lt0,ye,ent0,prec,ltout,2,keyerrt);
#if DEBUG
    if(*keyerrt==667) {
      fprintf(stderr,"This is worse. Bisection failed!\n");
      abort();
    }      
#endif


    return;
  }



} // namespace nuc_eos
