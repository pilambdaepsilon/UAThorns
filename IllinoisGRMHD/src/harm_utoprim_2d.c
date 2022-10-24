/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic 
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to 
    solve the relativistic magnetohydrodynamic equations of motion on a 
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model. 

    You are morally obligated to cite the following two papers in his/her 
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003, 
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

   
    Further, we strongly encourage you to obtain the latest version of 
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

utoprim_2d.c: 
---------------

    Uses the 2D method: 
       -- solves for two independent variables (W,v^2) via a 2D
          Newton-Raphson method 
       -- can be used (in principle) with a general equation of state. 

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want 
      to change this aspect of the code so that it still calculates the 
      velocity and so that you can floor the densities.  If you want to 
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();

******************************************************************************/

static const int NEWT_DIM=2;

// Declarations: 
static CCTK_REAL vsq_calc(CCTK_REAL W,CCTK_REAL &Bsq,CCTK_REAL &QdotBsq,CCTK_REAL &Qtsq,CCTK_REAL &Qdotn,CCTK_REAL &D);
static int Utoprim_new_body(CCTK_REAL U[], CCTK_REAL gcov[NDIM][NDIM], CCTK_REAL gcon[NDIM][NDIM], CCTK_REAL gdet,  CCTK_REAL prim[],long &n_iter);
static int general_newton_raphson( CCTK_REAL x[], int n, long &n_iter, void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [][NEWT_DIM], CCTK_REAL *, CCTK_REAL *, int,CCTK_REAL &,CCTK_REAL &,CCTK_REAL &,CCTK_REAL &,CCTK_REAL &),CCTK_REAL &Bsq,CCTK_REAL &QdotBsq,CCTK_REAL &Qtsq,CCTK_REAL &Qdotn,CCTK_REAL &D);
static void func_vsq( CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n,CCTK_REAL &Bsq,CCTK_REAL &QdotBsq,CCTK_REAL &Qtsq,CCTK_REAL &Qdotn,CCTK_REAL &D);
static CCTK_REAL x1_of_x0(CCTK_REAL x0, CCTK_REAL &Bsq, CCTK_REAL &QdotBsq, CCTK_REAL &Qtsq, CCTK_REAL &Qdotn, CCTK_REAL &D ) ;
static CCTK_REAL pressure_W_vsq(CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL &D) ;
static CCTK_REAL dpdW_calc_vsq(CCTK_REAL W, CCTK_REAL vsq);
static CCTK_REAL dpdvsq_calc(CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL &D);

/**********************************************************************/
/******************************************************************

  Utoprim_2d():
  
  -- Driver for new prim. var. solver.  The driver just translates
     between the two sets of definitions for U and P.  The user may 
     wish to alter the translation as they see fit.  Note that Greek
     indices run 0,1,2,3 and Latin indices run 1,2,3 (spatial only).


              /  rho u^t            \
         U =  |  T^t_t + rho u^t    |  sqrt(-det(g_{\mu\nu}))
              |  T^t_i              |
              \  B^i                /

             /    rho        \
     P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /


   Arguments:
       U[NPR]    = conserved variables (current values on input/output);
       gcov[NDIM][NDIM] = covariant form of the metric ;
       gcon[NDIM][NDIM] = contravariant form of the metric ;
       gdet             = sqrt( - determinant of the metric) ;
       prim[NPR] = primitive variables (guess on input, calculated values on 
                                        output if there are no problems);
  
   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set 
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.  ;

******************************************************************/

int Utoprim_2d(CCTK_REAL U[NPR], CCTK_REAL gcov[NDIM][NDIM], CCTK_REAL gcon[NDIM][NDIM], 
               CCTK_REAL gdet, CCTK_REAL prim[NPR], long &n_iter)
{

  CCTK_REAL U_tmp[NPR], prim_tmp[NPR];
  int i, ret; 
  CCTK_REAL alpha;

  if( U[0] <= 0. ) { 
    return(-100);
  }

  /* First update the primitive B-fields */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] / gdet ;

  /* Set the geometry variables: */
  alpha = 1.0/sqrt(-gcon[0][0]);
  
  /* Transform the CONSERVED variables into the new system */
  U_tmp[RHO] = alpha * U[RHO] / gdet;
  U_tmp[YECON]  = alpha * (U[YECON])  / gdet ;	//for Ye
  U_tmp[UU]  = alpha * (U[UU] - U[RHO])  / gdet ;
  for( i = UTCON1; i <= UTCON3; i++ ) {
    U_tmp[i] = alpha * U[i] / gdet ;
  }
  for( i = BCON1; i <= BCON3; i++ ) {
    U_tmp[i] = alpha * U[i] / gdet ;
  }

  /* Transform the PRIMITIVE variables into the new system */
  for( i = 0; i < BCON1; i++ ) {
    prim_tmp[i] = prim[i];
  }
  prim_tmp[YECON] = prim[YECON];		//for Ye
  for( i = BCON1; i <= BCON3; i++ ) {
    prim_tmp[i] = alpha*prim[i];
  }

  ret = Utoprim_new_body(U_tmp, gcov, gcon, gdet, prim_tmp,n_iter);

  /* Transform new primitive variables back if there was no problem : */ 
  if( ret == 0 || ret == 5 || ret==101 ) {
    for( i = 0; i < BCON1; i++ ) {
      prim[i] = prim_tmp[i];
    }
    prim[YECON] = prim_tmp[YECON];	//note were not doing it for B-fields here, but it will be done in Uprim_newbody
  }

  return( ret ) ;

}


/**********************************************************************/
/**********************************************************************************

  Utoprim_new_body():

     -- Attempt an inversion from U to prim using the initial guess prim.

     -- This is the main routine that calculates auxiliary quantities for the 
        Newton-Raphson routine. 

  -- assumes that 
             /  rho gamma        \
         U = |  alpha T^t_\mu  |
             \  alpha B^i        /



               /    rho        \
    prim = |    uu         |
               | \tilde{u}^i   |
               \  alpha B^i   /


return:  (i*100 + j)  where 
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used) 
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the 
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence 
                   (occurrence of "nan" or "+/-inf" ;
         
         j = 0 -> success 
             1 -> failure: some sort of failure in Newton-Raphson; 
             2 -> failure: utsq<0 w/ initial p[] guess;
         3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1 
             5 -> failure: rho,uu <= 0 ;

**********************************************************************************/

static int Utoprim_new_body(CCTK_REAL U[NPR], CCTK_REAL gcov[NDIM][NDIM], 
                            CCTK_REAL gcon[NDIM][NDIM], CCTK_REAL gdet,  CCTK_REAL prim[NPR], long &n_iter)
{

  CCTK_REAL x_2d[NEWT_DIM];
  CCTK_REAL QdotB,Bcon[NDIM],Bcov[NDIM],Qcov[NDIM],Qcon[NDIM],ncov[NDIM],ncon[NDIM],Qsq,Qtcon[NDIM];
  CCTK_REAL rho0,y_e,u,p,w,gammasq,gamma,gtmp,W_last,W,utsq,vsq;
  int i,j, n, retval, i_increase;

  n = NEWT_DIM ;

  // Assume ok initially:
  retval = 0;

  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  Bcon[0] = 0. ;
  for(i=1;i<4;i++) Bcon[i] = U[BCON1+i-1] ;

  lower_g(Bcon,gcov,Bcov) ;

  for(i=0;i<4;i++) Qcov[i] = U[QCOV0+i] ;
  raise_g(Qcov,gcon,Qcon) ;


  CCTK_REAL Bsq = 0. ;
  for(i=1;i<4;i++) Bsq += Bcon[i]*Bcov[i] ;

  QdotB = 0. ;
  for(i=0;i<4;i++) QdotB += Qcov[i]*Bcon[i] ;
  CCTK_REAL QdotBsq = QdotB*QdotB ;

  ncov_calc(gcon,ncov) ;
  raise_g(ncov,gcon,ncon);

  CCTK_REAL Qdotn = Qcon[0]*ncov[0] ;

  Qsq = 0. ;
  for(i=0;i<4;i++) Qsq += Qcov[i]*Qcon[i] ;

  CCTK_REAL Qtsq = Qsq + Qdotn*Qdotn ;

  CCTK_REAL D = U[RHO] ;
  CCTK_REAL Ye_t = U[YECON];			//Yet

  /* calculate W from last timestep and use for guess */
  utsq = 0. ;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++) utsq += gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1] ;


  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) { 
    utsq = fabs(utsq);
  }
  if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval) ;
  }

  gammasq = 1. + utsq ;
  gamma  = sqrt(gammasq);
    
  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) . 
  rho0 = D / gamma ;
  u = prim[UU] ;
  p = pressure_rho0_u(rho0,u) ;
  w = rho0 + u + p ;

  W_last = w*gammasq ;


  // Make sure that W is large enough so that v^2 < 1 : 
  i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*Bsq ) 
            - QdotBsq*(2.*W_last + Bsq) ) <= W_last*W_last*(Qtsq-Bsq*Bsq))
         && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }
  
  // Calculate W and vsq: 
  x_2d[0] =  fabs( W_last );
  x_2d[1] = x1_of_x0( W_last , Bsq,QdotBsq,Qtsq,Qdotn,D) ;
  retval = general_newton_raphson( x_2d, n, n_iter, func_vsq, Bsq,QdotBsq,Qtsq,Qdotn,D) ;  

  W = x_2d[0];
  vsq = x_2d[1];
    
  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (W == FAIL_VAL) ) {
    retval = retval*100+1;
    return(retval);
  }
  else{
    if(W <= 0. || W > W_TOO_BIG) {
      retval = 3;
      return(retval) ;
    }
  }

  // Calculate v^2:
  if( vsq >= 1. ) {
    vsq = 1.-2.e-16;
    //retval = 4;
    //return(retval) ;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  gtmp = sqrt(1. - vsq);
  gamma = 1./gtmp ;
  rho0 = D * gtmp;
  y_e = Ye_t/D;

  w = W * (1. - vsq) ;
  p = pressure_rho0_w(rho0,w) ;
  u = w - (rho0 + p) ; // u = rho0 eps, w = rho0 h

  if( (rho0 <= 0.) || (u <= 0.) ) {
    // User may want to handle this case differently, e.g. do NOT return upon 
    // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:

    retval = 5;
    //return(retval) ;
  }

  /*
    if(retval==5 && fabs(u)<1e-16) {
    u = fabs(u);
    CCTK_VInfo(CCTK_THORNSTRING,"%e\t%e\t%e",1.0-w/(rho0 + p),rho0,p);
    retval=0;
    }
  */

  prim[RHO] = rho0 ;
  prim[UU] = u ;
  prim[YECON] = y_e ;

  for(i=1;i<4;i++)  Qtcon[i] = Qcon[i] + ncon[i] * Qdotn;
  for(i=1;i<4;i++) prim[UTCON1+i-1] = gamma/(W+Bsq) * ( Qtcon[i] + QdotB*Bcon[i]/W ) ;
    
  /* set field components */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;


  /* done! */
  return(retval) ;

}


/**********************************************************************/ 
/****************************************************************************
   vsq_calc(): 
    
      -- evaluate v^2 (spatial, normalized velocity) from 
            W = \gamma^2 w 

****************************************************************************/
static CCTK_REAL vsq_calc(CCTK_REAL W,CCTK_REAL &Bsq,CCTK_REAL &QdotBsq,CCTK_REAL &Qtsq,CCTK_REAL &Qdotn,CCTK_REAL &D)
{
  CCTK_REAL Wsq,Xsq;
    
  Wsq = W*W ;
  Xsq = (Bsq + W) * (Bsq + W);

  return(  ( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq) );
}


/********************************************************************

  x1_of_x0(): 
           
    -- calculates v^2 from W  with some physical bounds checking;
    -- asumes x0 is already physical
    -- makes v^2 physical  if not;

*********************************************************************/

static CCTK_REAL x1_of_x0(CCTK_REAL x0, CCTK_REAL &Bsq, CCTK_REAL &QdotBsq, CCTK_REAL &Qtsq, CCTK_REAL &Qdotn, CCTK_REAL &D ) 
{
  CCTK_REAL vsq;
  CCTK_REAL dv = 1.e-15;
  
  vsq = fabs(vsq_calc(x0,Bsq,QdotBsq,Qtsq,Qdotn,D)) ; // guaranteed to be positive 


  return( ( vsq > 1. ) ? (1.0 - dv) : vsq   ); 

}

/********************************************************************

  validate_x(): 
           
    -- makes sure that x[0,1] have physical values, based upon 
       their definitions:
    
*********************************************************************/

static void validate_x(CCTK_REAL x[2], CCTK_REAL x0[2] ) 
{
  
  CCTK_REAL dv = 1.e-15;

  /* Always take the absolute value of x[0] and check to see if it's too big:  */ 
  x[0] = fabs(x[0]);
  x[0] = (x[0] > W_TOO_BIG) ?  x0[0] : x[0];
  

  x[1] = (x[1] < 0.) ?   0.       : x[1];  /* if it's too small */
  x[1] = (x[1] > 1.) ?  (1. - dv) : x[1];  /* if it's too big   */

  return;

}

/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static int general_newton_raphson( CCTK_REAL x[], int n, long &n_iter,
                                   void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], 
                                                  CCTK_REAL [][NEWT_DIM], CCTK_REAL *, 
                                                  CCTK_REAL *, int,CCTK_REAL &,CCTK_REAL &,CCTK_REAL &,CCTK_REAL &,CCTK_REAL &),CCTK_REAL &Bsq,CCTK_REAL &QdotBsq,CCTK_REAL &Qtsq,CCTK_REAL &Qdotn,CCTK_REAL &D)
{
  CCTK_REAL f, df, dx[NEWT_DIM], x_old[NEWT_DIM];
  CCTK_REAL resid[NEWT_DIM], jac[NEWT_DIM][NEWT_DIM];
  CCTK_REAL errx, x_orig[NEWT_DIM];
  int    id, i_extra, doing_extra;

  int   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;

  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, n, Bsq,QdotBsq,Qtsq,Qdotn,D);  /* returns with new dx, f, df */
      

    /* Save old values before calculating the new: */
    errx = 0.;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    /* Make the newton step: */
    for( id = 0; id < n ; id++) {
      x[id] += dx[id]  ;
    }

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    validate_x( x, x_old ) ;


    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) 
        || (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)

  /*  Check for bad untrapped divergences : */
  if( (std::isfinite(f)==0) ||  (std::isfinite(df)==0) ) {
    return(2);
  }


  if( fabs(errx) > MIN_NEWT_TOL){
    //CCTK_VInfo(CCTK_THORNSTRING,"%d %e %e %e %e",n_iter,f,df,errx,MIN_NEWT_TOL);
    return(1);
  } 
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }

  return(0);

}



/**********************************************************************/
/*********************************************************************************
   func_vsq(): 

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=W,vsq here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
*********************************************************************************/

static void func_vsq(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[], 
                     CCTK_REAL jac[][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n,
                     CCTK_REAL &Bsq,CCTK_REAL &QdotBsq,CCTK_REAL &Qtsq,CCTK_REAL &Qdotn,CCTK_REAL &D)
{

  
  CCTK_REAL  W, vsq, Wsq, p_tmp, dPdvsq, dPdW;
  CCTK_REAL t11;
  CCTK_REAL t16;
  CCTK_REAL t18;
  CCTK_REAL t2;
  CCTK_REAL t21;
  CCTK_REAL t23;
  CCTK_REAL t24;
  CCTK_REAL t25;
  CCTK_REAL t3;
  CCTK_REAL t35;
  CCTK_REAL t36;
  CCTK_REAL t4;
  CCTK_REAL t40;
  CCTK_REAL t9;

  // vv TESTING vv
  //  CCTK_REAL D,gtmp,gamma,rho0,w,p,u;
  // ^^ TESTING ^^

  W = x[0];
  vsq = x[1];
  
  Wsq = W*W;

  // vv TESTING vv
  /*
    D = U[RHO] ;
    gtmp = sqrt(1. - vsq);
    gamma = 1./gtmp ;
    rho0 = D * gtmp;

    w = W * (1. - vsq) ;
    p = pressure_rho0_w(rho0,w) ;
    u = w - (rho0 + p) ;

    if(u<=0 && 1==1) {
    vsq = 0.9999999 * (1.0-(rho0+p)/W);

    w = W * (1. - vsq) ;
    p = pressure_rho0_w(rho0,w) ;
    u = w - (rho0 + p) ;

    //CCTK_VInfo(CCTK_THORNSTRING,"%e check",u);
    }
  */
  // ^^ TESTING ^^

  
  p_tmp  = pressure_W_vsq( W, vsq , D);
  dPdW   = dpdW_calc_vsq( W, vsq );
  dPdvsq = dpdvsq_calc( W, vsq, D );

  // These expressions were calculated using Mathematica, but made into efficient 
  // code using Maple.  Since we know the analytic form of the equations, we can 
  // explicitly calculate the Newton-Raphson step: 

  t2 = -0.5*Bsq+dPdvsq;
  t3 = Bsq+W;
  t4 = t3*t3;
  t9 = 1/Wsq;
  t11 = Qtsq-vsq*t4+QdotBsq*(Bsq+2.0*W)*t9;
  t16 = QdotBsq*t9;
  t18 = -Qdotn-0.5*Bsq*(1.0+vsq)+0.5*t16-W+p_tmp;
  t21 = 1/t3;
  t23 = 1/W;
  t24 = t16*t23;
  t25 = -1.0+dPdW-t24;
  t35 = t25*t3+(Bsq-2.0*dPdvsq)*(QdotBsq+vsq*Wsq*W)*t9*t23;
  t36 = 1/t35;
  dx[0] = -(t2*t11+t4*t18)*t21*t36;
  t40 = (vsq+t24)*t3;
  dx[1] = -(-t25*t11-2.0*t40*t18)*t21*t36;
  //detJ = t3*t35; // <- set but not used...
  jac[0][0] = -2.0*t40;
  jac[0][1] = -t4;
  jac[1][0] = t25;
  jac[1][1] = t2;
  resid[0] = t11;
  resid[1] = t18;

  

  *df = -resid[0]*resid[0] - resid[1]*resid[1];

  *f = -0.5 * ( *df );

}


/********************************************************************** 
 ********************************************************************** 
   
 The following routines specify the equation of state.  All routines 
  above here should be indpendent of EOS.  If the user wishes 
  to use another equation of state, the below functions must be replaced 
  by equivalent routines based upon the new EOS. 

  **********************************************************************
  **********************************************************************/

/**********************************************************************/
/********************************************************************** 
  pressure_W_vsq():  
  
        -- Gamma-law equation of state;
        -- pressure as a function of W, vsq, and D:
**********************************************************************/
static CCTK_REAL pressure_W_vsq(CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL &D) 
{
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL gtmp;
  gtmp = 1. - vsq;
  
  return(  (gamma_th /* <- Should be local polytropic Gamma factor */  - 1.) * ( W * gtmp  -  D * sqrt(gtmp) ) / gamma_th /* <- Should be local polytropic Gamma factor */   );

}


/**********************************************************************/
/********************************************************************** 
  dpdW_calc_vsq(): 
 
      -- partial derivative of pressure with respect to W;
**********************************************************************/
static CCTK_REAL dpdW_calc_vsq(CCTK_REAL W, CCTK_REAL vsq)
{
  DECLARE_CCTK_PARAMETERS;
  return( (gamma_th /* <- Should be local polytropic Gamma factor */  - 1.) * (1. - vsq) /  gamma_th /* <- Should be local polytropic Gamma factor */  ) ;

}

/**********************************************************************/
/********************************************************************** 
  dpdvsq_calc(): 
 
      -- partial derivative of pressure with respect to vsq
**********************************************************************/
static CCTK_REAL dpdvsq_calc(CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL &D)
{
  DECLARE_CCTK_PARAMETERS;
  return( (gamma_th /* <- Should be local polytropic Gamma factor */  - 1.) * ( 0.5 * D / sqrt(1.-vsq)  - W  ) / gamma_th /* <- Should be local polytropic Gamma factor */   ) ;
}


/****************************************************************************** 
             END   OF   UTOPRIM_2D.C
******************************************************************************/



