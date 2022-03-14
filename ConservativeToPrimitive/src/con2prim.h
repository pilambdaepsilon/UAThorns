#ifndef con2prim
#define con2prim

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "con2prim_config.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define NDIM 4
#define NPR  9

#define RHO 0
#define v1_cov 1
#define v2_cov 2
#define v3_cov 3
#define EPS 4

#define v1_contra 1
#define v2_contra 2
#define v3_contra 3

#define D 0
#define S1_cov 1
#define S2_cov 2
#define S3_cov 3
#define TAU 4

#define B1_con 5
#define B2_con 6
#define B3_con 7

#define YE 8
#define TEMP 9
#define PRESS 10
#define WLORENTZ 11
#define ENT 12
#define A_BAR 13
#define MU_HAT 14

#define DEBUG_NUNO 0
#define ONE_OVER_SQRT_4PI 0.282094791773878143474039725780
#define MAX_NTAB 16001
#define IMAX(a,b) ( a>b ? a : b ) 
#define IMIN(a,b) ( a<b ? a : b ) 
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define DBL_EPSILON_C2P 1e-15

static const int C2P_PHI=0,C2P_PSI=1,C2P_GXX=2,C2P_GXY=3,C2P_GXZ=4,C2P_GYY=5,C2P_GYZ=6,C2P_GZZ=7,
  C2P_LAPM1=8,C2P_SHIFTX=9,C2P_SHIFTY=10,C2P_SHIFTZ=11,C2P_GUPXX=12,C2P_GUPYY=13,C2P_GUPZZ=14;
static const int C2P_GUPXY=15,C2P_GUPXZ=16,C2P_GUPYZ=17;
static const int C2P_LAPSE=0,C2P_PSI2=1,C2P_PSI4=2,C2P_PSI6=3,C2P_PSIM4=4,C2P_LAPSEINV=5;
static const int C2P_RHOB=0,C2P_PRESSURE=1,C2P_VX=2,C2P_VY=3,C2P_VZ=4,
  C2P_BX_CENTER=5,C2P_BY_CENTER=6,C2P_BZ_CENTER=7,C2P_BX_STAGGER=8,C2P_BY_STAGGER=9,C2P_BZ_STAGGER=10,
  C2P_VXR=11,C2P_VYR=12,C2P_VZR=13,C2P_VXL=14,C2P_VYL=15,C2P_VZL=16,C2P_YE=17,C2P_EPS=18,C2P_TEMP=19;  //<-- Be _sure_ to define MAXNUMVARS appropriately!
static const int C2P_RHOSTAR=0,C2P_STILDEX=1,C2P_STILDEY=2,C2P_STILDEZ=3,C2P_TAUENERGY=4,C2P_YET=5;
//static const int C2P_RHO=0, C2P_UU=1,C2P_UTCON1=2, UTCON2=3, UTCON3=4, BCON1=5, BCON2=6, BCON3=7, YECON=8;
static const int UU=1,UTCON1=2,UTCON2=3,UTCON3=4,BCON1=5,BCON2=6,BCON3=7,YECON=8;
static const int QCOV0=1,QCOV1=2,QCOV2=3,QCOV3=4;
static const CCTK_REAL UTSQ_TOO_BIG =1.e20;
static const CCTK_REAL C2P_NEWT_TOL    =1.0e-10;
static const CCTK_REAL C2P_MIN_NEWT_TOL    =1.0e-10;
static const int C2P_MAX_NEWT_ITER=30;
static const int C2P_EXTRA_NEWT_ITER=0;

static const CCTK_REAL W_TOO_BIG    =1.e20;


  // constants, variables and structures needed
  // by standalone version

  // some constants (required by EOS_Omni
  // in standalone version)
  extern const double kBerg;
  extern const double amu;

  // conversion factors between cgs and M_Sun = c = G = 1
  // (temperature between K and MeV)
  // see EOS_Omni/doc/units.py
  extern const double rho_gf;
  extern const double press_gf;
  extern const double eps_gf;
  extern const double temp_gf;

  // Inverses of the numbers above
  extern const double inv_rho_gf;
  extern const double inv_press_gf;
  extern const double inv_eps_gf;
  extern const double inv_temp_gf;


  // EOS related quantities

  // nuc eos

//  struct nuc_eos_vars {
//    int nrho;
//    int ntemp;
//    int nye;
//    double *fourtables;
//    double *logrho;
//    double *logtemp;
//    double *yes;
//    double energy_shift;
//    double dtemp, dtempi;
//    double drho, drhoi;
//    double dye, dyei;
//    // min and max values
//    double eos_rhomax, eos_rhomin;
//    double eos_tempmin, eos_tempmax;
//    double eos_yemin, eos_yemax;
//    double eos_epsmin, eos_epsmax;
//    int *ivs_short;
//  };



struct metric {
    double lo[4][4];
    double up[4][4];
    double lo_det;
    double up_det;
    double lo_sqrt_det;
};


struct c2p_steer {
  CCTK_INT c2p_method;
  CCTK_INT c2p_method_backup;
  bool use_c2p_method_backup;
  CCTK_INT eoskey;
  CCTK_INT eoskey_polytrope;
  CCTK_INT key_temp;
  CCTK_REAL eos_prec;
  CCTK_REAL eos_rho_min;
  CCTK_REAL eos_rho_max;
  CCTK_REAL eos_temp_min;
  CCTK_REAL eos_temp_max;
  CCTK_REAL eos_eps_min;
  CCTK_REAL eos_eps_max;
  CCTK_REAL eos_press_min;
  CCTK_REAL eos_press_max;
  CCTK_REAL eos_ye_min;
  CCTK_REAL eos_ye_max;
  CCTK_INT evolve_T;
  CCTK_INT evolve_Ye;
  CCTK_REAL tol_x;
  CCTK_REAL tol_x_retry;
  CCTK_INT max_iterations;
  CCTK_INT extra_iterations;
  CCTK_REAL rho_atmo;
  CCTK_REAL rho_atmo_tol;
  CCTK_REAL T_atmo;
  CCTK_REAL Ye_atmo;
  CCTK_REAL P_atmo;
  CCTK_REAL eps_atmo;
  bool retain_B_atmo;
  CCTK_INT numprims;
  CCTK_INT numcons;
  bool enforce_v2;
  bool force_cold_beta_equil;
};

struct c2p_steer c2p;

struct c2p_report {
  bool failed;
  bool adjust_cons;
  char err_msg[200];
  int count;
  bool retry;
  int c2p_keyerr;
  int nEOScalls;
};

/* This is visible everywhere */
extern struct c2p_steer c2p;
extern struct c2p_report extern_report;

// functions



void con2prim_MHD_(double * prim, const double * con, CCTK_REAL g_up[4][4],
      CCTK_REAL g_lo[4][4], const bool key_excise, const bool c2p_grace,
      struct c2p_report * rep, const int cac_ite, const int index, const double xco, const double yco, const double zco, double * METRIC, double * METRIC_PHYS, double * METRIC_LAP_PSI4);

void NR_3D(struct c2p_report * c2p_rep, const double S_squared, const double BdotS,
  const double B_squared, const double * con, double * prim, const double g_con[4][4],
  const double g_cov[4][4], const double tol_x, int SAFEGUESS, int stepsize);

void palenzuela(struct c2p_report * c2p_rep, const double S_squared, const double BdotS,
      const double B_squared, const double * con, double * prim, const double g_con[4][4],
      const double g_cov[4][4], const double tol_x, bool use_epsmin, const int index,double xco, double yco, double zco,int iteration);

void NR_2D_Noble(struct c2p_report * c2p_rep, const double S_squared, const double BdotS,
      const double B_squared, const double * con, double * prim, const double g_con[4][4],
      const double g_cov[4][4], const double tol_x);

void HARM(struct c2p_report * c2p_rep, const double S_squared, const double BdotS,
      const double B_squared, const double * con, double * prim, CCTK_REAL g_con[4][4],
      CCTK_REAL g_cov[4][4], const double tol_x,double xco, double yco, double zco,const int iteration, const int index, double * METRIC, double * METRIC_PHYS, double * METRIC_LAP_PSI4);

double zbrent(double (*func)(double, double *, struct c2p_report *, bool, double *), double * param, double * temp_guess, double x1, double x2, double tol, struct c2p_report * c2p_rep, bool use_epsmin);

void newman(struct c2p_report * c2p_rep, const double S_squared,
          const double BdotS, const double B_squared, const double * con, double * prim,
	      const double g_cov[4][4], const double g_con[4][4], const double tol_x, double xco, double yco, double zco,int iteration, double * METRIC_LAP_PSI4);

void NR_2D_rhoT(struct c2p_report * c2p_rep, const double S_squared,
      const double BdotS, const double B_squared, const double * con, double * prim,
    const double glo[4][4], const double tol_x, int SAFEGUESS, int stepsize);

void NR_2D_WT(struct c2p_report * c2p_rep, const double S_squared,
      const double BdotS, const double B_squared, const double * con, double * prim,
    const double glo[4][4], const double tol_x, int SAFEGUESS, int stepsize);

void EOS_press_ent_abar(int eoskey, int keytemp, double rho, double * eps,
        double * temp, double ye, double * prs, double * ent, double * abar, int * keyerr);


void EOS_EP_dEdr_dEdt_dPdr_dPdt(double * x, const double * con, double * Eprim,
    double * Pprim, double * dEdrho, double * dEdt, double * dPdrho, double * dPdt, int stepsize);

void EOS_EP_dEdr_dEdt_dPdr_dPdt_2D(const double rho2D, const double temp2D, const double * con, double * Eprim,
    double * Pprim, double * dEdrho, double * dEdt, double * dPdrho, double * dPdt, int stepsize);

void EP_dEdW_dEdZ_dEdT(double * Eprim, double * Pprim, double * dEdW, double * dEdZ, double * dEdT,
    double * dpdrho, double * dpdT, double * x, const double * con, int stepsize);

void EP_dPdW_dPdZ_dPdT(double * Eprim, double * Pprim, double * dPdW, double * dPdZ, double * dPdT,
    double * dpdrho, double * dpdT, double * x, const double * con, int stepsize);

void EOS_P_from_hrho_dPdrho_dPdeps(const double rho, const double enth, const double * con,
          double * temp_guess, double * eps, double * press, double * dPdrho,
          double * dPdeps, double * entr, double * abar);

#ifdef __cplusplus
}
#endif

#endif
