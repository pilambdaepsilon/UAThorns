#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "consts.h"
#include "rnsid_util.h"
#include "equil_util.h"

#define ARRAY_END 1
#define FREE_ARG char*

/*************************************************************************
* Allocate memory for matrix and 3d arrays (tensor)
*************************************************************************/
void AllocationError(char error_text[])
{
	fprintf(stderr,"RNS RUNTIME ERROR: %s\n Exiting the system ",error_text);
	exit(1);
}


double **array_allocate(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+ARRAY_END)*sizeof(double*)));
	if (!m) AllocationError("allocation failure 1 in matrix()");
	m += ARRAY_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+ARRAY_END)*sizeof(double)));
	if (!m[nrl]) AllocationError("allocation failure 2 in matrix()");
	m[nrl] += ARRAY_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


double ***tensor_allocate(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+ARRAY_END)*sizeof(double**)));
	if (!t) AllocationError("allocation failure 1 in tensor_allocate()");
	t += ARRAY_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+ARRAY_END)*sizeof(double*)));
	if (!t[nrl]) AllocationError("allocation failure 2 in tensor_allocate()");
	t[nrl] += ARRAY_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+ARRAY_END)*sizeof(double)));
	if (!t[nrl][ncl]) AllocationError("allocation failure 3 in tensor_allocate()");
	t[nrl][ncl] += ARRAY_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}



/*************************************************************************
* Free memory of double Tensor and ARRAY 
*************************************************************************/

void array_free(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by array_allocate() */
{
	free((FREE_ARG) (m[nrl]+ncl-ARRAY_END));
	free((FREE_ARG) (m+nrl-ARRAY_END));
}


void tensor_free(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a double tensor_allocate allocated by tensor_allocate() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-ARRAY_END));
	free((FREE_ARG) (t[nrl]+ncl-ARRAY_END));
	free((FREE_ARG) (t+nrl-ARRAY_END));
}



/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate using four points around xb.       */  
/*************************************************************************/
double interp_4(double xp[5], 
                double yp[5], 
                int    np ,
                double xb)
{ 
 int k=1;      /* index of 1st point */ 
 double y;     /* intermediate value */


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



/*************************************************************************
* Interpolation between two different grids
*************************************************************************/
void grid_interp(double **old, 
                 double *s_gp, 
                 double *mu, 
                 double r_e, 
                 int nx,
                 int ny,
                 int nz,  
                 double *x_grid,
                 double *y_grid,
                 double *z_grid,
                 int i,
                 int j,
                 int k, 
                 double *new,
                 int sign) {

int s,
    m,
    s_nearest,
    m_nearest,                /* nearest points in interpolation */
    k_s,                      /* first s point in interpolation */
    k_m;                      /* first s point in interpolation */


double r_c,                   /* r of cartesian x,y,z point */
       s_c,                   /* s of cartesian x,y,z point */
       mu_c,                  /* mu of cartesian x,y,z point */
       s_4[5],                /* s of the 4 nearest points */
       mu_4[5],               /* mu of the 4 nearest points */
       old_s[5],              /* old at 4 nearest constant s points */
       old_m[5];              /* old at 4 nearest constant mu points */  


	    
      r_c = sqrt(   SQ(x_grid[i-1+nx*(j-1+ny*(k-1))])   
                       + SQ(y_grid[i-1+nx*(j-1+ny*(k-1))]) 
                       + SQ(z_grid[i-1+nx*(j-1+ny*(k-1))]) ); 
      s_c = r_c/(r_e+r_c);
      
      if(r_c==0.0) 
        mu_c = 0.0;
      else
        mu_c = fabs(z_grid[i-1+nx*(j-1+ny*(k-1))])/r_c;
 
      s_nearest = 0; m_nearest = 0;

      hunt(s_gp, SDIV, s_c, &s_nearest);
      hunt(mu, MDIV, mu_c, &m_nearest); 
  
      k_s = IMIN(IMAX((s_nearest)-(4-1)/2,1),SDIV+1-4);
      k_m = IMIN(IMAX((m_nearest)-(4-1)/2,1),MDIV+1-4);
      	       
      for(s=1;s<=4;s++) 
            s_4[s] = s_gp[k_s-1+s];

      for(m=1;m<=4;m++) 
            mu_4[m] = mu[k_m-1+m];
 
      for(s=1;s<=4;s++) {
          for(m=1;m<=4;m++) {
                  old_s[m] = old[k_s-1+s][k_m-1+m];
          }
          old_m[s] = interp_4(mu_4, old_s, 4, mu_c);  
      }

     if (z_grid[i-1+nx*(j-1+ny*(k-1))]<0.0)                               
              (*new) = (1.0*sign)*interp_4(s_4, old_m, 4, s_c); 
     else 
              (*new) = interp_4(s_4, old_m, 4, s_c); 
}



/*************************************************************************
* Interpolation between two different grids - all
*************************************************************************/
void grid_interp_all( double *s_gp, 
                      double *mu, 
                      double  r_e, 
                      int     nx,
                      int     ny,
                      int     nz, 
                      double *x_grid,
                      double *y_grid,
                      double *z_grid, 
                      int i,
                      int j,
                      int k, 
                      double **nu, 
                      double **B, 
                      double **alpha, 
                      double **omega, 
                      double **nu_dr, 
                      double **B_dr, 
                      double **alpha_dr, 
                      double **omega_dr, 
                      double **nu_dth, 
                      double **B_dth, 
                      double **alpha_dth, 
                      double **omega_dth,
                      double **rho_0, 
                      double **energy, 
                      double **pressure, 
                      double *nu_c, 
                      double *B_c, 
                      double *alpha_c, 
                      double *omega_c, 
                      double *nu_dr_c, 
                      double *B_dr_c, 
                      double *alpha_dr_c, 
                      double *omega_dr_c,                        
                      double *nu_dth_c, 
                      double *B_dth_c, 
                      double *alpha_dth_c, 
                      double *omega_dth_c,
                      double *rho_0_c, 
                      double *energy_c, 
                      double *pressure_c,
                      double *distance_c,
		      double **Omega_diff,
		      double *Omega_diff_c)
{  
  
      grid_interp( nu, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, nu_c, 1);
      
      grid_interp( B, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, B_c,  1);  

      grid_interp( alpha, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, alpha_c,  1);

      grid_interp( omega, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, omega_c,  1);

      grid_interp( nu_dr, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, nu_dr_c, 1);

      grid_interp( B_dr, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, B_dr_c, 1);

      grid_interp( alpha_dr, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, alpha_dr_c, 1);

      grid_interp( omega_dr, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, omega_dr_c, 1);
		   
      grid_interp( nu_dth, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, nu_dth_c, -1);
      
      grid_interp( B_dth, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, B_dth_c, -1);

      grid_interp( alpha_dth, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, alpha_dth_c, -1);

      grid_interp( omega_dth, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, omega_dth_c, -1);
		   	   
      grid_interp( rho_0, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, rho_0_c, 1);
      
      grid_interp( energy, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, energy_c, 1);

      grid_interp( pressure, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, pressure_c, 1);       	       

      grid_interp( Omega_diff, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, Omega_diff_c, 1);

}

