//Routines for cold beta equilibrium

void load_beta_equilibrium_file( const char beta_equil_file[],
               double log_rho0_table[MAX_NTAB],
               double Y_e_table[MAX_NTAB],
               int *n_tab_beta)
{
 int i;                    /* counter */

 /*constants to convert from cgs to cactus units c=G=M_sun=1.0 */
 CCTK_REAL const cactusM= (5.028916268544129e-34);    /*  1/g  */
 CCTK_REAL const cactusL= (6.772400341316594e-06);    /*  1/cm */
 CCTK_REAL const cactusT= (2.0303145448833407e5);    /*  1/s  */
 CCTK_REAL const cactusV= (1.0/(cactusL*cactusL*cactusL));
 double rho0,               /* density */
        ye;                /* electron fraction */

 FILE *f_beta;              /* pointer to beta_equil_file */


    /* OPEN FILE TO READ */

    if((f_beta=fopen(beta_equil_file,"r")) == NULL ) {
       CCTK_VERROR("cannot open beta-equil. file:  %s\n",beta_equil_file);
    }

    /* READ NUMBER OF TABULATED POINTS */

    fscanf(f_beta,"%d\n",n_tab_beta);

    /* READ ENERGY DENSITY, P, H, N0 AND CONVERT TO CACTUS UNITS */

    for(i=1;i<=(*n_tab_beta);i++) {
       fscanf(f_beta,"%lf %lf \n",&rho0,&ye) ;
       log_rho0_table[i]=log10(rho0*cactusM*cactusV);     /* multiply by C^2 to get energy density */
       if(ye <= 0.036){
           Y_e_table[i] = 0.036;
       }
       else{
           Y_e_table[i]=ye;
       }
    }
}

/***************************************************************************/
/* Routine that locates nearest grid point for a given value.              */
/***************************************************************************/
void hunt_c2p(double xx[], int n, double x, int *jlo)
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

/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate using four points around xb.       */  
/*************************************************************************/
double interp_c2p(double xp[], 
              double yp[], 
              int    np ,
              double xb, 
              int    *n_nearest_pt)
{ 
 int k,        /* index of 1st point */
     m=4;      /* degree of interpolation */ 
 
 double y;     /* intermediate value */

 hunt_c2p(xp,np,xb,n_nearest_pt);

 k=IMIN(IMAX((*n_nearest_pt)-(m-1)/2,1),np+1-m);

 if( xb==xp[k] ||  xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) 
    xb += DBL_EPSILON_C2P;

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


