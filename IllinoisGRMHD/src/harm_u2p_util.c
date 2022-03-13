
/*
  -------------------------------------------------------------------------------
  Copyright 2005 Scott C. Noble, Charles F. Gammie, 
  Jonathan C. McKinney, and Luca Del Zanna


  This file is part of PVS-GRMHD.

  PVS-GRMHD is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  PVS-GRMHD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with PVS-GRMHD; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  -------------------------------------------------------------------------------
*/

// Function prototypes for this file:
void primtoU_g( CCTK_REAL prim[], CCTK_REAL gcov[][NDIM], CCTK_REAL gcon[][NDIM], CCTK_REAL gdet,  CCTK_REAL U[] );
static void ucon_calc_g(CCTK_REAL prim[],CCTK_REAL gcov[][NDIM],CCTK_REAL gcon[][NDIM],CCTK_REAL ucon[]);
static void raise_g(CCTK_REAL vcov[], CCTK_REAL gcon[][NDIM], CCTK_REAL vcon[]);
static void lower_g(CCTK_REAL vcon[], CCTK_REAL gcov[][NDIM], CCTK_REAL vcov[]);
static void ncov_calc(CCTK_REAL gcon[][NDIM],CCTK_REAL ncov[]) ;
static void bcon_calc_g(CCTK_REAL prim[],CCTK_REAL ucon[],CCTK_REAL ucov[],CCTK_REAL ncov[],CCTK_REAL bcon[]); 
static CCTK_REAL pressure_rho0_u(CCTK_REAL rho0, CCTK_REAL u);
static CCTK_REAL pressure_rho0_w(CCTK_REAL rho0, CCTK_REAL w);
int gamma_calc_g(CCTK_REAL *pr, CCTK_REAL gcov[NDIM][NDIM], CCTK_REAL *gamma);

/********************************************************************** 
   primtoU_g(): 

       -- calculates the conserved variables from the primitive variables 
            and the metric;
       -- assumes that the conserved and primitive variables are defined ala HARM:

              /  rho u^t           \
         U =  |  T^t_\mu + rho u^t |  sqrt(-det(g_{\mu\nu}))
              |   B^i              |
	      \   Ye_tilde         /

             /    rho        \
     P =     |   uu          |
             | \tilde{u}^i   |
             |   B^i         |
	     \   Ye	     /

**************************************************************************/


void primtoU_g(
               CCTK_REAL prim[NPR],       /* primitive variables */
               CCTK_REAL gcov[NDIM][NDIM],    /* covariant (index dn) form of metric */
               CCTK_REAL gcon[NDIM][NDIM],    /* contravariant (index up) form of metric */
               CCTK_REAL gdet,                /* sqrt of -1 times det(g_{\mu \nu}) */
               CCTK_REAL U[NPR]           /* matrix of derivatives */
               ) {
  int i ;
  CCTK_REAL rho0 ;
  CCTK_REAL y_e ;
  static CCTK_REAL ucon[NDIM],ucov[NDIM],bcon[NDIM],bcov[NDIM],ncov[NDIM] ;
  CCTK_REAL gamma,n_dot_b,bsq,u,p,w, alpha ;

    
  /* Calculate auxiliary quantities: */
  alpha = 1.0/sqrt(-gcon[0][0]);

  ucon_calc_g(prim,gcov,gcon,ucon) ;
  lower_g(ucon,gcov,ucov) ;
  ncov_calc(gcon,ncov) ;

  gamma = -ncov[0]*ucon[0] ;

  bcon_calc_g(prim,ucon,ucov,ncov,bcon) ;
  lower_g(bcon,gcov,bcov) ;

  n_dot_b = 0. ;
  for(i=0;i<4;i++) n_dot_b += ncov[i]*bcon[i] ;
  bsq = 0. ;
  for(i=0;i<4;i++) bsq += bcov[i]*bcon[i] ;

  rho0 = prim[RHO] ;
  y_e = prim[YECON] ;
  u = prim[UU] ;
  p = pressure_rho0_u(rho0,u) ;
  w = rho0 + u + p ;

  // Now set the conserved variables themselves, using HARM's definition:
  U[RHO] = ucon[0]*rho0 ;
  U[YECON] = U[RHO]*y_e ;

  for( i = 0; i < 4; i++) {
    U[QCOV0+i] = gamma*(w + bsq)*ucov[i] 
      - (p + bsq/2.)*ncov[i] 
      + n_dot_b*bcov[i] ;

    U[QCOV0+i] /= alpha;
  }

  U[QCOV0] = U[QCOV0] + U[RHO];
  U[BCON1] = prim[BCON1] ;
  U[BCON2] = prim[BCON2] ;
  U[BCON3] = prim[BCON3] ;

  for(i = 0; i < NPR; i++ ) {
    U[i] *= gdet;
  }

  return ;
}

/********************************************************************** 
  ucon_calc_g(): 
    
       -- calculates the contravariant (up) components of the four-velocity
          given the primitive variables, of which the velocity is 
          \tilde{u}^i = \gamma v^j  where v^j is the velocity of the 
          flow w.r.t a normal observer to the coordinates;

       -- also requires the metric and inverse metric;

       -- assumes:

             /    rho        \
     P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /

******************************************************************/
static void ucon_calc_g(CCTK_REAL prim[NPR],CCTK_REAL gcov[NDIM][NDIM],CCTK_REAL gcon[NDIM][NDIM],
                 CCTK_REAL ucon[NDIM])
{
  CCTK_REAL u_tilde_con[4] ;
  CCTK_REAL u_tilde_sq ;
  CCTK_REAL gamma,lapse ;
  int i,j;
    
  u_tilde_con[0] = 0. ;
  u_tilde_con[1] = prim[UTCON1] ;
  u_tilde_con[2] = prim[UTCON2] ;
  u_tilde_con[3] = prim[UTCON3] ;

  u_tilde_sq = 0. ;
  for(i=0;i<NDIM;i++)
    for(j=0;j<NDIM;j++)
      u_tilde_sq += gcov[i][j]*u_tilde_con[i]*u_tilde_con[j] ;
  u_tilde_sq = fabs(u_tilde_sq) ;

  gamma = sqrt(1. + u_tilde_sq) ;

  lapse = sqrt(-1./gcon[0][0]) ;

  for(i=0;i<NDIM;i++) ucon[i] = u_tilde_con[i] - lapse*gamma*gcon[0][i] ;

  return ;
}

/********************************************************************** 
    raise_g():
 
         -- calculates the contravariant form of a covariant tensor, 
            using the inverse of the metric;
******************************************************************/
static void raise_g(CCTK_REAL vcov[NDIM], CCTK_REAL gcon[NDIM][NDIM], CCTK_REAL vcon[NDIM])
{
  int i,j;

  for(i=0;i<NDIM;i++) {
    vcon[i] = 0. ;
    for(j=0;j<NDIM;j++) 
      vcon[i] += gcon[i][j]*vcov[j] ;
  }

  return ;
}

/********************************************************************** 
     lower_g():
  
          -- calculates the ocvariant form of a contravariant tensor 
             using the metric;
******************************************************************/
static void lower_g(CCTK_REAL vcon[NDIM], CCTK_REAL gcov[NDIM][NDIM], CCTK_REAL vcov[NDIM])
{
  int i,j;

  for(i=0;i<NDIM;i++) {
    vcov[i] = 0. ;
    for(j=0;j<NDIM;j++) 
      vcov[i] += gcov[i][j]*vcon[j] ;
  }

  return ;
}

/********************************************************************** 
     ncov_calc(): 

         -- calculates the covariant form of the normal vector to our 
            spacelike hypersurfaces ala the ADM formalism.

         -- requires the inverse metric;
******************************************************************/
static void ncov_calc(CCTK_REAL gcon[NDIM][NDIM],CCTK_REAL ncov[NDIM]) 
{
  CCTK_REAL lapse ;
  int i;

  lapse = sqrt(-1./gcon[0][0]) ;

  ncov[0] = -lapse ;
  for( i = 1; i < NDIM; i++) { 
    ncov[i] = 0. ;
  }

  return ;
}

/********************************************************************** 
    bcon_calc_g(): 
  
        -- using the primitive variables, contra-/co-variant 4-vel., 
           and covariant normal vector, calculate the contravariant 
           form of the magnetic 4-vector b^\mu (the small "b" in HARM);
       -- assumes:

             /    rho        \
     P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /
******************************************************************/
static void bcon_calc_g(CCTK_REAL prim[NPR],CCTK_REAL ucon[NDIM],CCTK_REAL ucov[NDIM],
                 CCTK_REAL ncov[NDIM],CCTK_REAL bcon[NDIM]) 
{
  static CCTK_REAL Bcon[NDIM] ;
  CCTK_REAL u_dot_B ;
  CCTK_REAL gamma ;
  int i ;

  // Bcon = \mathcal{B}^\mu  of the paper:
  Bcon[0] = 0. ;
  for(i=1;i<NDIM;i++) Bcon[i] = -ncov[0] * prim[BCON1+i-1] ;

  u_dot_B = 0. ;
  for(i=0;i<NDIM;i++) u_dot_B += ucov[i]*Bcon[i] ;

  gamma = -ucon[0]*ncov[0] ;
  for(i=0;i<NDIM;i++) bcon[i] = (Bcon[i] + ucon[i]*u_dot_B)/gamma ;
}


/********************************************************************** 
    gamma_calc_g(): 
  
        -- using the primitive variables, contra-/co-variant 4-vel., 
           and covariant normal vector, calculate the contravariant 
           form of the magnetic 4-vector b^\mu (the small "b" in HARM);

       -- assumes:

             /    rho        \
     P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /
******************************************************************/
int gamma_calc_g(CCTK_REAL *pr, CCTK_REAL gcov[NDIM][NDIM], CCTK_REAL *gamma)
{
  CCTK_REAL utsq ;

  utsq =    gcov[1][1]*pr[UTCON1]*pr[UTCON1]
    + gcov[2][2]*pr[UTCON2]*pr[UTCON2]
    + gcov[3][3]*pr[UTCON3]*pr[UTCON3]
    + 2.*(  gcov[1][2]*pr[UTCON1]*pr[UTCON2]
            + gcov[1][3]*pr[UTCON1]*pr[UTCON3]
            + gcov[2][3]*pr[UTCON2]*pr[UTCON3]) ;

  if(utsq<0.0){
    if(fabs(utsq)>1E-10){ // then assume not just machine precision
      return (1);
    }
    else utsq=1E-10; // set floor
  }

  *gamma = sqrt(1. + utsq) ;

  return(0) ;
}


/**************************************************
  The following functions assume a Gamma-law EOS:
***************************************************/

/* 
   pressure as a function of rho0 and u 
   this is used by primtoU and Utoprim_?D
*/
static CCTK_REAL pressure_rho0_u(CCTK_REAL rho0, CCTK_REAL u)
{
  DECLARE_CCTK_PARAMETERS;
  return((gamma_th /* <- Should be local polytropic Gamma factor */  - 1.)*u) ;
}


  
/* 
   pressure as a function of rho0 and w = rho0 + u + p 
   this is used by primtoU and Utoprim_1D
*/
static CCTK_REAL pressure_rho0_w(CCTK_REAL rho0, CCTK_REAL w)
{
  DECLARE_CCTK_PARAMETERS;
  return((gamma_th /* <- Should be local polytropic Gamma factor */ -1.)*(w - rho0)/gamma_th /* <- Should be local polytropic Gamma factor */ ) ;
}


