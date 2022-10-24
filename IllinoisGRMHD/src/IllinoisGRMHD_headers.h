// To safeguard against double-including this header file:
#ifndef ILLINOISGRMHD_HEADERS_H_
#define ILLINOISGRMHD_HEADERS_H_

#include "ConservativeToPrimitive.h"

#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define SQR(x) ((x) * (x))
#define ONE_OVER_SQRT_4PI 0.282094791773878143474039725780
#define DEBUG_NUNO 0
#define DEBUG_NUNO_PREPOS 0

#define VERR_DEF_PARAMS __LINE__, __FILE__, CCTK_THORNSTRING

// The order here MATTERS, as we assume that GUPXX+1=GUPYY, etc.
static const int PHI=0,PSI=1,GXX=2,GXY=3,GXZ=4,GYY=5,GYZ=6,GZZ=7,
  LAPM1=8,SHIFTX=9,SHIFTY=10,SHIFTZ=11,GUPXX=12,GUPYY=13,GUPZZ=14,
  NUMVARS_FOR_METRIC_FACEVALS=15; //<-- Be _sure_ to set this correctly, or you'll have memory access bugs!

// These are not used for facevals in the reconstruction step, but boy are they useful anyway.
static const int GUPXY=15,GUPXZ=16,GUPYZ=17,
  NUMVARS_FOR_METRIC=18; //<-- Be _sure_ to set this correctly, or you'll have memory access bugs!

// The order here MATTERS, and must be consistent with the order in the in_prims[] array in driver_evaluate_MHD_rhs.C.
static const int RHOB=0,PRESSURE=1,VX=2,VY=3,VZ=4,
  BX_CENTER=5,BY_CENTER=6,BZ_CENTER=7,BX_STAGGER=8,BY_STAGGER=9,BZ_STAGGER=10,
  VXR=11,VYR=12,VZR=13,VXL=14,VYL=15,VZL=16,YE=17,EPS=18,TEMP=19,WLORENTZ=20,ENTROPY=21,MAXNUMVARS=22;  //<-- Be _sure_ to define MAXNUMVARS appropriately!

static const int UT=0,UX=1,UY=2,UZ=3;

// The "I" suffix denotes interpolation. In other words, these
//    definitions are used for interpolation ONLY. The order here
//    matters as well!
static const int SHIFTXI=0,SHIFTYI=1,SHIFTZI=2,GUPXXI=3,GUPXYI=4,GUPXZI=5,GUPYYI=6,GUPYZI=7,GUPZZI=8,
  PSII=9,LAPM1I=10,A_XI=11,A_YI=12,A_ZI=13,LAPSE_PSI2I=14,LAPSE_OVER_PSI6I=15,MAXNUMINTERP=16;

// Again, the order here MATTERS, since we assume in the code that, e.g., smallb[0]=b^t, smallb[3]=b^z, etc.
static const int SMALLBT=0,SMALLBX=1,SMALLBY=2,SMALLBZ=3,SMALLB2=4,NUMVARS_SMALLB=5;

// Again, the order here MATTERS, since we assume in the code that, CONSERV[STILDEX+1] = \tilde{S}_y
static const int RHOSTAR=0,STILDEX=1,STILDEY=2,STILDEZ=3,TAUENERGY=4,YET=5,NUM_CONSERVS=6;

static const int LAPSE=0,PSI2=1,PSI4=2,PSI6=3,PSIM4=4,LAPSEINV=5,NUMVARS_METRIC_AUX=6;
#define SET_LAPSE_PSI4(array_name,METRIC)   {                   \
      array_name[LAPSE] = METRIC[LAPM1]+1.0;                    \
      array_name[PSI2]  = exp(2.0*METRIC[PHI]);                 \
      array_name[PSI4]  = SQR(array_name[PSI2]);                \
      array_name[PSI6]  = array_name[PSI4]*array_name[PSI2];    \
      array_name[PSIM4]  = 1.0/array_name[PSI4];                \
      array_name[LAPSEINV]  = 1.0/array_name[LAPSE];            \
  }

// Keeping track of ghostzones between routines is a nightmare, so
//   we instead attach ghostzone info to each gridfunction and set
//   the ghostzone information correctly within each routine.
struct gf_and_gz_struct {
  CCTK_REAL *gf;
  int gz_lo[4],gz_hi[4];
};

#define MAX_EOS_PARAMS 10
struct eos_struct {
  int neos;
  CCTK_REAL rho_tab[MAX_EOS_PARAMS],P_tab[MAX_EOS_PARAMS],eps_tab[MAX_EOS_PARAMS],k_tab[MAX_EOS_PARAMS],gamma_tab[MAX_EOS_PARAMS];
  CCTK_REAL gamma_th;
  CCTK_REAL K_poly;
};

struct output_stats {
  int font_fixed,vel_limited,failure_checker;
  long n_iter;
};


// FIXME: For cosmetic purposes, we might want to make everything either zero-offset or one-offset, instead of a mixture.
const int kronecker_delta[4][3] = { { 0,0,0 },
                                    { 1,0,0 },
                                    { 0,1,0 },
                                    { 0,0,1 } };

/* PUBLIC FUNCTIONS, USED OUTSIDE IllinoisGRMHD AS WELL */
void IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(const int eos_key, CCTK_REAL yemin, CCTK_REAL yemax, const int already_computed_physical_metric_and_inverse,CCTK_REAL *U,
		struct output_stats &stats,eos_struct &eos, CCTK_REAL *METRIC,CCTK_REAL g4dn[4][4],CCTK_REAL g4up[4][4], CCTK_REAL *TUPMUNU,CCTK_REAL *TDNMUNU,CCTK_REAL *CONSERVS, const int index);

void IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij
(const cGH *cctkGH,const int *cctk_lsh,
 CCTK_REAL *gxx,CCTK_REAL *gxy,CCTK_REAL *gxz,CCTK_REAL *gyy,CCTK_REAL *gyz,CCTK_REAL *gzz,CCTK_REAL *alp,
 CCTK_REAL *gtxx,CCTK_REAL *gtxy,CCTK_REAL *gtxz,CCTK_REAL *gtyy,CCTK_REAL *gtyz,CCTK_REAL *gtzz,
 CCTK_REAL *gtupxx,CCTK_REAL *gtupxy,CCTK_REAL *gtupxz,CCTK_REAL *gtupyy,CCTK_REAL *gtupyz,CCTK_REAL *gtupzz,
 CCTK_REAL *phi,CCTK_REAL *psi,CCTK_REAL *lapm1);

void IllinoisGRMHD_set_symmetry_gzs_staggered(const cGH *cctkGH, const int *cctk_lsh,CCTK_REAL *X,CCTK_REAL *Y,CCTK_REAL *Z,  CCTK_REAL *gridfunc,
                                              CCTK_REAL *gridfunc_syms,int stagger_x,int stagger_y,int stagger_z);

#endif // ILLINOISGRMHD_HEADERS_H
