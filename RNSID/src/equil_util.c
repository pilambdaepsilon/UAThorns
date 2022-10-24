/*****************************************************************************
equil_util.c	UTILITY PROGRAMS USED TO COMPUTE THE EQUILIBRIUM
       deriv_s() etc.: returns the derivative of a function wrt s, etc
*****************************************************/


#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <stdlib.h>
#include "consts.h"
#include "rnsid_util.h"
#include "equil_util.h"


/***************************************************************************/
/* Routine that locates nearest grid point for a given value.              */
/***************************************************************************/
void hunt(double xx[], int n, double x, int *jlo)
{ 
	int jm,jhi,inc,ascnd;

	ascnd=(xx[n] > xx[1]);
	if (*jlo <= 0 || *jlo > n) {
		*jlo=0;
		jhi=n+1;
	} else {
		inc=1;
		if (x >= xx[*jlo] == ascnd) {
			if (*jlo == n) return;
			jhi=(*jlo)+1;
			while (x >= xx[jhi] == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n) {
					jhi=n+1;
					break;
				}
			}
		} else {
			if (*jlo == 1) {
				*jlo=0;
				return;
			}
			jhi=(*jlo);
			*jlo -= 1;
			while (x < xx[*jlo] == ascnd) {
				jhi=(*jlo);
				inc += inc;
				*jlo=jhi-inc;
				if (*jlo < 1) {
					*jlo=0;
					break;
				}
			}
		}
	}
	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if (x > xx[jm] == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
}

/*C*/
/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate using four points around xb.       */  
/*************************************************************************/
double interp(double xp[], 
              double yp[], 
              int    np ,
              double xb, 
              int    *n_nearest_pt)
{ 
 int k,        /* index of 1st point */
     m=4;      /* degree of interpolation */ 
 
 double y;     /* intermediate value */

 hunt(xp,np,xb,n_nearest_pt);

 k=IMIN(IMAX((*n_nearest_pt)-(m-1)/2,1),np+1-m);

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


/*******************************************************************/
/* Returns the derivative w.r.t. s                                 */ 
/*******************************************************************/
double deriv_s(double **f,int s, int m)
{ 
 double d_temp;

 switch(s) { 
            case 1    : d_temp=(f[s+1][m]-f[s][m])/DS;
                        break;

            case SDIV   : d_temp=(f[s][m]-f[s-1][m])/DS;
                          break;
      
            default     : d_temp=(f[s+1][m]-f[s-1][m])/(2.0*DS);
                          break; 
 } 
 
 return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s                                 */ 
/*******************************************************************/
double deriv_ss(double **f,int s, int m)
{ 
 double d_temp;

 switch(s) { 
/*
            case 1    : d_temp=(-f[s+3][m]+4.0*f[s+2][m]-5.0*f[s+1][m]
                                                     +2.0*f[s][m] )/SQ(DS);
                        break;

            case SDIV   : d_temp=(2.0*f[s][m]-5.0*f[s-1][m]+4.0*f[s-2][m]
                                                            -f[s-3][m] )/SQ(DS);
                          break;
      
            default     : d_temp=(f[s+1][m]-2.0*f[s][m]+f[s-1][m])/SQ(DS);
                          break; 
*/
/* 
           case 1    : d_temp=(f[s][m]-2.0*f[s+1][m]+f[s+2][m])/(2.0*SQ(DS));
                       break;

           case 2    : d_temp=(2.0*f[s-1][m]-3.0*f[s][m]+f[s+2][m])/(4.0*SQ(DS));
                       break;

           case SDIV-1 : d_temp=(2.0*f[s+1][m]-3.0*f[s][m]+f[s-2][m])
                                                                 /(4.0*SQ(DS));
                         break;

           case SDIV   : d_temp=(f[s][m]-2.0*f[s-1][m]+f[s-2][m])/(2.0*SQ(DS));
                          break;
      
           default     : d_temp=(f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*SQ(DS));
                         break; 
*/
 
           case 1    : s=4;
                       d_temp=(f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*SQ(DS));
                       break;

           case 2    : s=4;
                       d_temp=(f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*SQ(DS));
                       break;

           case 3    : s=4;
                       d_temp=(f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*SQ(DS));
                       break;

           case SDIV-1 : s=SDIV-2;
                         d_temp=(f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*SQ(DS));
                         break;

           case SDIV   :  s=SDIV-2;
                          d_temp=(f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*SQ(DS));
                          break;
      
           default     : d_temp=(f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*SQ(DS));
                         break; 


 } 
 return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. mu                                */ 
/*******************************************************************/
double deriv_m(double **f,int s, int m)
{
 double d_temp;

 switch(m) { 
            case 1    : d_temp=(f[s][m+1]-f[s][m])/DM;
                        break; 

            case MDIV   : d_temp=(f[s][m]-f[s][m-1])/DM;
                          break;
      
            default     : d_temp=(f[s][m+1]-f[s][m-1])/(2.0*DM);
                          break; 
 } 
 return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s                                 */ 
/*******************************************************************/
double deriv_mm(double **f,int s, int m)
{ 
 double d_temp;

 switch(m) { 
            case 1    : m=2;
                        d_temp=(f[s][m+1]-2.0*f[s][m]+f[s][m-1])/SQ(DM);
                        break;

            case MDIV   : m=MDIV-1;
                          d_temp=(f[s][m+1]-2.0*f[s][m]+f[s][m-1])/SQ(DM);
                          break;

            default     : d_temp=(f[s][m+1]-2.0*f[s][m]+f[s][m-1])/SQ(DM);
                          break; 
 } 
 return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s and mu                          */ 
/*******************************************************************/
double deriv_sm(double **f,int s, int m)
{
 double d_temp;

 switch(s) {
     case 1 : if(m==1) {   
               d_temp=(f[s+1][m+1]-f[s][m+1]-f[s+1][m]+f[s][m])/(DM*DS);
              }else{
                if(m==MDIV) {
                 d_temp=(f[s+1][m]-f[s][m]-f[s+1][m-1]+f[s][m-1])/(DM*DS);
                }else{         
                   d_temp=(f[s+1][m+1]-f[s+1][m-1]-f[s][m+1]+f[s][m-1])/
                                                                (2.0*DM*DS);
                }
              }
              break;

     case SDIV : if(m==1) {   
               d_temp=(f[s][m+1]-f[s][m]-f[s-1][m+1]+f[s-1][m])/(DM*DS);
              }else{
                if(m==MDIV) {
                 d_temp=(f[s][m]-f[s-1][m]-f[s][m-1]+f[s-1][m-1])/(DM*DS);
                }else{         
                   d_temp=(f[s][m+1]-f[s][m-1]-f[s-1][m+1]+f[s-1][m-1])/
                                                                (2.0*DM*DS);
                }
             }
             break;
  
     default : if(m==1) {   
               d_temp=(f[s+1][m+1]-f[s-1][m+1]-f[s+1][m]+f[s-1][m])/(2.0*DM*DS);
              }else{
                if(m==MDIV) {
                 d_temp=(f[s+1][m]-f[s-1][m]-f[s+1][m-1]+f[s-1][m-1])/
                                                                (2.0*DM*DS);
                }else{         
                  d_temp=(f[s+1][m+1]-f[s-1][m+1]-f[s+1][m-1]+f[s-1][m-1])/
                                                                (4.0*DM*DS);
                }
             }
             break;
     }

  return d_temp;

}




// ********************************************************************
// Compute the Legendre polinomial of order (n) at (X) using recursion
//
//  P_0(x) = 1,   P_1(x) = x  
//
//  P_n(x) = (2*n-1)/n * P_{n-1}(x) - (n-1)/n P_{n-2} (x)
//
// ********************************************************************
 
double legendre(int n, double  x) {
  int i;           /* counter */

  double p,        /* Legendre polynomial of order n */
         p_1,      /*    "         "      "    "   n-1*/
         p_2;      /*    "         "      "    "   n-2 */

  p_2=1.0;
  p_1=x;

 if(n>=2) { 
  for(i=2;i<=n;i++){
     p=(x*(2.0*i-1.0)*p_1 - (i-1.0)*p_2)/i;
     p_2=p_1;
     p_1=p;
  }
  return p;
 } else { 
    if (n==1) return p_1;
      else return p_2;
   }
}

//  double p,p1,p2 ; // Legendre polynomial of order n,n-1,n.2  
//  double di;
//  switch (n) {
//  case 0:
//    return 1.0;
//  case 1:
//    return x;
//    } 
//  p2 = 0;
//  p1 = x;
//  for (int i=2;i<=n;i++) {
//    di = (double)(i);
//    p  = 1/di*( (2.0*di-1.0)*p1 - (di-1.0)* p2   );
//    p2 = p1;
//    p1 = p ;
//    }
//  return p;
//
//}
//
//
// ************************************************************************
// Compute the Associate Legendre polinomial of order P_l^m(x) 
// using recursion. Assume m>=0,m<=l, -1<=x<=1
//
//  P_0^0(x)     = 1.0
//  P_{m}^{m}(x) = - (2*m-1) sqrt(1-x**2) P_{m-1}^{m-1}(x)  
//
//  P_{m+1}^m(x) = x * (2m+1) P_m^m(x)  
//
//  P_l^m(x) = 1/(l-m) ( (2*l-1) * x P_{l-1}^m(x) - (l-1+m) * P_{l-2}^m(x)  
//
//*************************************************************************

double plgndr(int l, int m, double x)
{
	double fact,pll,pmm,pmmp1,somx2;
	int i,ll;

	if (m < 0 || m > l || fabs(x) > 1.0)
		printf("Bad arguments in routine PLGNDR");
	pmm=1.0;
	if (m > 0) {
		somx2=sqrt((1.0-x)*(1.0+x));
		fact=1.0;
		for (i=1;i<=m;i++) {
			pmm *= -fact*somx2;
			fact += 2.0;
		}
	}
	if (l == m)
		return pmm;
	else {
		pmmp1=x*(2*m+1)*pmm;
		if (l == (m+1))
			return pmmp1;
		else {
			for (ll=(m+2);ll<=l;ll++) {
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm=pmmp1;
				pmmp1=pll;
			}
			return pll;
		}
	}
}




/*******************************************************************/
/* (based on rtsec) */
/*******************************************************************/
double rtsec_G(double (*func)(double, double, double), 
               double Gamma_P, 
               double x1, 
               double x2, 
               double xacc, 
               double ee,
	       double eos_k)

{
 int j;
 double fl,f,dx,swap, xl,rts;
 
 fl=(*func)(x1,Gamma_P,eos_k)-ee;
 f=(*func)(x2,Gamma_P,eos_k)-ee;

 if(fabs(fl)<fabs(f)) {
   rts=x1;
   xl=x2;
   swap=fl;
   fl=f;
   f=swap;
 } else {
         xl=x1;
         rts=x2;
        }

 // #define MAXIT 1000 
 for(j=1;j<=MAXIT;j++) {
    dx=(xl-rts)*f/(f-fl);
    xl=rts;
    fl=f;
    rts += dx;
    f=(*func)(rts,Gamma_P,eos_k)-ee;

    if(fabs(dx)<xacc||f==0.0) return rts;
  }
 
 printf("Maximum number of iterations exceeded in rtsec");  
 return 0.0;
}


/**************************************************************************/
/* (based on zbrent) */
/**************************************************************************/
double zbrent_diff(double (*func)(double, double, double, double, double, 
                                  double, double, double), 
                   double r_e,
                   double rho_equator_h,
                   double gama_equator_h,
                   double omega_equator_h,
                   double rho_pole_h,
                   double gama_pole_h,  
                   double x1, 
                   double x2, 
                   double tol,
		   double A_diff)
{
	int    iter;

	double a=x1,
               b=x2,
               c=x2,
               d,
               e,
               min1,
               min2;

	double fa=(*func)(a, r_e, rho_equator_h, gama_equator_h, 
                          omega_equator_h, rho_pole_h, gama_pole_h, A_diff),
               fb=(*func)(b, r_e, rho_equator_h, gama_equator_h, 
                          omega_equator_h, rho_pole_h, gama_pole_h, A_diff),
               fc,
               p,
               q, 
               r,
               s,
               tol1,
	       xm;
	       
	double return_value;

	if ((fa > 0.0 && fb>0.0) || (fa<0.0 && fb<0.0)) {
                printf("Root must be bracketed in ZBRENT_DIFF\n");
				/* exit(1); */
	}
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc>0.0) || (fb<0.0 && fc<0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) 
		  { return_value=b;
    		    return return_value;
		  }  

		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0)  q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1, xm);
		fb=(*func)(b, r_e, rho_equator_h, gama_equator_h,  
                           omega_equator_h, rho_pole_h, gama_pole_h, A_diff);
	} 
	printf("Maximum number of iterations exceeded in ZBRENT_DIFF\n");
	exit(1);

	return return_value;
}


/**************************************************************************/
/* (based on zbrent) */
/**************************************************************************/
double zbrent_rot(double (*func)(double, double, double, double, double, 
                                  double, double, double), 
                   double r_e,
                   double rhogp, 
                   double omegagp, 
                   double sgp, 
                   double mugp, 
                   double Omega_c,
                   double x1, 
                   double x2, 
                   double tol,
		   double A_diff)
{
	int    iter;

	double a=x1,
               b=x2,
               c=x2,
               d,
               e,
               min1,
               min2;

	double fa=(*func)(a, r_e, rhogp, omegagp, 
                          sgp, mugp, Omega_c, A_diff),
               fb=(*func)(b, r_e, rhogp, omegagp, 
                          sgp, mugp, Omega_c, A_diff),
               fc,
               p,
               q, 
               r,
               s,
               tol1,
               xm;

	double return_value;

	if ((fa > 0.0 && fb>0.0) || (fa<0.0 && fb<0.0)) {
                printf("Root must be bracketed in ZBRENT_ROT at s_gp=%4.3e, mu=%4.3e\n", sgp,mugp);
				exit(1);
	}
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc>0.0) || (fb<0.0 && fc<0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) 
		  { return_value= b;
		    return return_value;
		  }
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0)  q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1, xm);
                fb=(*func)(b, r_e, rhogp, omegagp, sgp, mugp, Omega_c, A_diff);
	} 
	printf("Maximum number of iterations exceeded in ZBRENT_ROT\n");
    exit(1);
	
	return return_value;

}


/**************************************************************************/
double diff_rotation(double x, 
                     double r_e, 
                     double rho_equator_h, 
                     double gama_equator_h, 
                     double omega_equator_h, 
                     double rho_pole_h, 
                     double gama_pole_h,
		     double A_diff)
{
 return (SQ(r_e)*(gama_equator_h+rho_equator_h-gama_pole_h-rho_pole_h)
        +log(1.0-SQ( (x-omega_equator_h)*exp(-SQ(r_e)*rho_equator_h) )))*
        SQ( 1.0 - SQ((x-omega_equator_h)*exp(-SQ(r_e)*rho_equator_h)))
        -SQ( (1.0/A_diff)*(x-omega_equator_h)*exp(-2.0*SQ(r_e)*rho_equator_h));
} 

/**************************************************************************/
double rotation_law( double x, 
                     double r_e, 
                     double rhogp, 
                     double omegagp, 
                     double sgp, 
                     double mugp, 
                     double Omega_c,
		     double A_diff)
{

 return (Omega_c -x)*( SQ(1.0-sgp) - SQ( sgp*(x-omegagp)*exp(-SQ(r_e)*rhogp) )
        *(1.0-mugp*mugp) ) - SQ((1.0/A_diff))*(x-omegagp)*sgp*sgp
        *(1.0-mugp*mugp)*exp(-2.0*SQ(r_e)*rhogp);
}
