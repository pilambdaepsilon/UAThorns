/* LOW RESOLUTION (USE THIS ON LOW-MEMORY MACHINES) */
/*
#define MDIV 81
#define SDIV 151
*/
/* OTHER POSSIBLE CHOICES */ 
/*
#define MDIV 161
#define SDIV 301
*/
/*
#define MDIV 121
#define SDIV 341
*/
/* 
#define MDIV 221
#define SDIV 341
*/
/*
#define MDIV 221
#define SDIV 541
*/
/*
#define MDIV 161
#define SDIV 841
*/
/*#define MDIV 161*/ /*201*/             /* Number of angular grid points */
/*#define SDIV 181*/                     /* Number of radial grid points */
/*#define SDIV 501*/ /* 1001 */
/*#define SDIV 401 */

/* MAXIMUM RESOLUTION (NO POINT IN GOING BEYOND THAT DUE TO GIBBS ERRORS) */
#define MDIV 301 
#define SDIV 601
#define DM (1.0/(MDIV-1.0))          /* spacing in mu direction */ 
#define RDIV 900                     /* grid point in RK integration */ 
#define SMAX 0.9999                  /* maximum value of s-coordinate */  
#define DS (SMAX/(SDIV-1.0))         /* spacing in s-direction */

/* ================================================================================
  Constant set from the value present in LORENE 
  Lorene/C++/Include/unites.h

  const double g_si = 6.6738E-11 ;     ///< Newton gravitational constant [SI]
  const double c_si = 2.99792458E+8 ;  ///< Velocity of light [m/s]
  const double kB_si = 1.3806488E-23 ; ///< Boltzmann constant [J/K]
  const double rhonuc_si = 1.66E+17 ;  ///< Nuclear density [kg/m3] (arbitrary)
  const double km_si = 1.E+3 ;         ///< One kilometer [m]
  const double msol_si = 1.9885E+30 ;  ///< Solar mass [kg]
  const double mev_si = 1.602176565E-13 ;   ///< One MeV [J]

  const double mu_si = 1.2566370614359173e-6 ;///<Magnetic vacuum permeability

  The Review of Particle Physics (2016)
  C. Patrignani et al. (Particle Data Group), Chin. Phys. C, 40, 100001 (2016). 
  http://www-pdg.lbl.gov/2016/reviews/rpp2016-rev-phys-constants.pdf
  http://www-pdg.lbl.gov/2016/reviews/rpp2016-rev-astrophysical-constants.pdf
  ==============================================================================  */

#define C 2.99792458e10              /* speed of light in vacuum */
#define G 6.6738e-8                  /* gravitational constant */ 
#define KAPPA (1.0e-15*C*C/G)        /* scaling factor */
#define KSCALE (KAPPA*G/(C*C*C*C))   /* another scaling factor */
#define MSUN 1.9885e33               /* Mass of Sun */
#define SQ(x) ((x)*(x))              /* square macro */
#define MB 1.66e-24                  /* baryon mass (conventional) */
#define RMIN 1.0e-15                 /* use approximate TOV equations when
					computing spherical star and r<RMIN */

/*constants to convert from cgs to cactus units c=G=M_sun=1.0 */
#define cactusM (5.028916268544129e-34)    /*  1/g  */
#define cactusL (6.772400341316594e-06)    /*  1/cm */
#define cactusT (2.0303145448833407e5)     /*  1/s  */
#define cactusV (1.0/(cactusL*cactusL*cactusL))

#define MAXIT 128                       /* MAX ITERATION IN SECANT */ 
#ifndef PI
#define PI (3.141592653589793)          /* what else */
#endif
#include <float.h>
#include <limits.h>

#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-15
#endif

#define MAX_NTAB 16001         /* MAX ELEMENTS IN EOS TAB                     */
#define MAX_ITERATION 2000     /* max number of iteration in iterate (euil.c) */
#define ITMAX 2000             /* MAX ITERATION IN BRENT  (equil_util.c)      */ 
#define EPS DBL_EPSILON
