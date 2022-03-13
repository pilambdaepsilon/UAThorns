/*-----------------------------------------------------------------------------
 *
 * Apply "Font Fix", of Font et al., ONLY when the primitives (con2prim) solver
 *     fails. Note that this version is OPTIMIZED for gamma=2 polytrope EOS only.
 *     Application of this fix GUARANTEES (on an analytic level, at least) 
 *     that con2prim will succeed.
 *
 * FIXME: There exist better fixes/strategies now. See McKinney et al's work.
 *
 -----------------------------------------------------------------------------*/

inline int font_fix_gamma_equals2(CCTK_REAL &u_x, CCTK_REAL &u_y, CCTK_REAL &u_z,CCTK_REAL *CONSERVS,CCTK_REAL *PRIMS,CCTK_REAL *METRIC_PHYS,CCTK_REAL *METRIC_LAP_PSI4,eos_struct &eos) {
  CCTK_REAL rhob;
  CCTK_REAL kpoly2 = 2.0*eos.K_poly;
  CCTK_REAL tol = 1.e-15;
  CCTK_REAL Bxbar = PRIMS[BX_CENTER]*ONE_OVER_SQRT_4PI;
  CCTK_REAL Bybar = PRIMS[BY_CENTER]*ONE_OVER_SQRT_4PI;
  CCTK_REAL Bzbar = PRIMS[BZ_CENTER]*ONE_OVER_SQRT_4PI;
  CCTK_REAL Bbar_x = METRIC_PHYS[GXX]*Bxbar + METRIC_PHYS[GXY]*Bybar + METRIC_PHYS[GXZ]*Bzbar;
  CCTK_REAL Bbar_y = METRIC_PHYS[GXY]*Bxbar + METRIC_PHYS[GYY]*Bybar + METRIC_PHYS[GYZ]*Bzbar;
  CCTK_REAL Bbar_z = METRIC_PHYS[GXZ]*Bxbar + METRIC_PHYS[GYZ]*Bybar + METRIC_PHYS[GZZ]*Bzbar;
  CCTK_REAL B2bar = Bxbar*Bbar_x + Bybar*Bbar_y + Bzbar*Bbar_z;
  CCTK_REAL Bbar = sqrt(B2bar);
  CCTK_REAL check_B_small = fabs(Bxbar)+fabs(Bybar)+fabs(Bzbar);
  if (check_B_small>0 && check_B_small<1.e-150) {
    // need to compute B2bar specially to prevent floating-point underflow
    CCTK_REAL Bmax = fabs(Bxbar);
    if (Bmax < fabs(Bybar)) Bmax=fabs(Bybar);
    if (Bmax < fabs(Bzbar)) Bmax=fabs(Bzbar);
    CCTK_REAL Bxtmp=Bxbar/Bmax, Bytemp=Bybar/Bmax, Bztemp=Bzbar/Bmax;
    CCTK_REAL B_xtemp=Bbar_x/Bmax, B_ytemp=Bbar_y/Bmax, B_ztemp=Bbar_z/Bmax;
    Bbar = sqrt(Bxtmp*B_xtemp + Bytemp*B_ytemp + Bztemp*B_ztemp)*Bmax;
  }

  CCTK_REAL BbardotS = Bxbar*CONSERVS[STILDEX] + Bybar*CONSERVS[STILDEY] + Bzbar*CONSERVS[STILDEZ];
  CCTK_REAL BbardotS2 = BbardotS*BbardotS;
  CCTK_REAL hatBbardotS = BbardotS/Bbar;
  if (Bbar<1.e-300) hatBbardotS = 0.0;
  CCTK_REAL Psim6 = 1.0/METRIC_LAP_PSI4[PSI6];

  // Limit hatBbardotS
  //CCTK_REAL max_gammav = 100.0;
  //CCTK_REAL rhob_max = CONSERVS[RHOSTAR]*Psim6;
  //CCTK_REAL hmax = 1.0 + kpoly2*rhob_max;
  //CCTK_REAL abs_hatBbardotS_max = sqrt(SQR(max_gammav)-1.0)*CONSERVS[RHOSTAR]*hmax;
  //if (fabs(hatBbardotS) > abs_hatBbardotS_max) {
  //   CCTK_REAL fac_reduce = abs_hatBbardotS_max/fabs(hatBbardotS);
  //   CCTK_REAL hatBbardotS_max = hatBbardotS*fac_reduce;
  //   CCTK_REAL Bbar_inv = 1.0/Bbar;
  //   CCTK_REAL hat_Bbar_x = Bbar_x*Bbar_inv;
  //   CCTK_REAL hat_Bbar_y = Bbar_y*Bbar_inv;
  //   CCTK_REAL hat_Bbar_z = Bbar_z*Bbar_inv;
  //   CCTK_REAL sub_fact = hatBbardotS_max - hatBbardotS;
  //   CONSERVS[STILDEX] += sub_fact*hat_Bbar_x;
  //   CONSERVS[STILDEY] += sub_fact*hat_Bbar_y;
  //   CONSERVS[STILDEZ] += sub_fact*hat_Bbar_z;
  //   hatBbardotS = hatBbardotS_max;
  //   BbardotS *= fac_reduce;
  //   BbardotS2 = BbardotS*BbardotS;
  //}

  CCTK_REAL sdots = METRIC_PHYS[GUPXX]*SQR(CONSERVS[STILDEX]) + METRIC_PHYS[GUPYY]*SQR(CONSERVS[STILDEY]) + METRIC_PHYS[GUPZZ]*SQR(CONSERVS[STILDEZ])
    + 2.0*( METRIC_PHYS[GUPXY]*CONSERVS[STILDEX]*CONSERVS[STILDEY] + METRIC_PHYS[GUPXZ]*CONSERVS[STILDEX]*CONSERVS[STILDEZ]
            + METRIC_PHYS[GUPYZ]*CONSERVS[STILDEY]*CONSERVS[STILDEZ]);

  if (sdots<1.e-300) {
    rhob = CONSERVS[RHOSTAR]*Psim6;
    u_x=0.0; u_y=0.0; u_z=0.0;
    return 0;
  }

  /*
    if (fabs(BbardotS2 - sdots*B2bar) > 1e-8) {
    CCTK_VInfo(CCTK_THORNSTRING,"(Bbar dot S)^2, Bbar^2 * sdotS, %e %e",SQR(BbardotS),sdots*B2bar);
    CCTK_VInfo(CCTK_THORNSTRING,"Cauchy-Schwartz inequality is violated!");
    }
  */

  // Initial guess for W, S_fluid and rhob
  CCTK_REAL W0 = sqrt( SQR(hatBbardotS) + SQR(CONSERVS[RHOSTAR]) ) * Psim6;
  CCTK_REAL Sf20 = (SQR(W0)*sdots + BbardotS2*(B2bar + 2.0*W0))/SQR(W0+B2bar);
  CCTK_REAL rhob0 = CONSERVS[RHOSTAR]*Psim6/sqrt(1.0+Sf20/SQR(CONSERVS[RHOSTAR]));
  CCTK_REAL W=W0,Sf2=Sf20,rhob1=rhob0;

  //****************************************************************
  //                          FONT FIX
  // Impose Font fix when HARM primitives solver fails to find
  //   acceptable set of primitives.
  //****************************************************************
  bool fontcheck=true;

  int itcount = 0, maxits=300;
  while(fontcheck && itcount < maxits) {
    itcount++;
    W0 = W;
    Sf20 = Sf2;
    rhob0 = rhob1;

    // first, find rhob for the given S_fluid^2
    CCTK_REAL h = 1.0 + kpoly2*rhob0;
    rhob1 = CONSERVS[RHOSTAR]*Psim6/sqrt(1.0+Sf20/SQR(CONSERVS[RHOSTAR]*h));
    while( fabs(rhob1-rhob0) > rhob1*tol) {
      rhob0 = rhob1;
      h = 1.0 + kpoly2*rhob0;
      rhob1 = CONSERVS[RHOSTAR]*Psim6/sqrt(1.0+Sf20/SQR(CONSERVS[RHOSTAR]*h));
    }

    h = 1.0 + kpoly2*rhob1;
    W = sqrt( Sf20 + SQR(CONSERVS[RHOSTAR]*h))*Psim6;
    Sf2 = (SQR(W)*sdots + BbardotS2*(B2bar + 2.0*W))/SQR(W+B2bar);
    if ( fabs(W-W0) < W*tol && fabs(Sf20-Sf2) < Sf2*tol) fontcheck=false;
  }
  
  if (itcount>=maxits) {
    // Increase tol and try again 
    maxits*=100;
    tol *=10.0;
    itcount = 0;
    fontcheck=true;
    while(fontcheck && itcount < maxits) {
      itcount++;
      W0 = W;
      Sf20 = Sf2;
      rhob0 = rhob1;

      // first, find rhob for the given S_fluid^2
      CCTK_REAL h = 1.0 + kpoly2*rhob0; 
      rhob1 = CONSERVS[RHOSTAR]*Psim6/sqrt(1.0+Sf20/SQR(CONSERVS[RHOSTAR]*h));
      while( fabs(rhob1-rhob0) > rhob1*tol) {
        rhob0 = rhob1;
        h = 1.0 + kpoly2*rhob0;
        rhob1 = CONSERVS[RHOSTAR]*Psim6/sqrt(1.0+Sf20/SQR(CONSERVS[RHOSTAR]*h));
      }

      h = 1.0 + kpoly2*rhob1;
      W = sqrt( Sf20 + SQR(CONSERVS[RHOSTAR]*h))*Psim6;
      Sf2 = (SQR(W)*sdots + BbardotS2*(B2bar + 2.0*W))/SQR(W+B2bar);
      if ( fabs(W-W0) < W*tol && fabs(Sf20-Sf2) < Sf2*tol) fontcheck=false;
    }
  }
  //************************************************************************************************************** 

  if(fontcheck==true) {
    return 1;
  }

  // Font fix works, now compute u_i
  rhob = rhob1;
  CCTK_REAL h = 1.0+kpoly2*rhob;
  CCTK_REAL gammav = CONSERVS[RHOSTAR]*Psim6/rhob;
  CCTK_REAL rhosh = CONSERVS[RHOSTAR]*h;
  CCTK_REAL fac1 = METRIC_LAP_PSI4[PSI6]*BbardotS/(gammav*rhosh);
  CCTK_REAL fac2 = 1.0/(rhosh + METRIC_LAP_PSI4[PSI6]*B2bar/gammav);
  u_x = fac2*(CONSERVS[STILDEX] + fac1*Bbar_x);
  u_y = fac2*(CONSERVS[STILDEY] + fac1*Bbar_y);
  u_z = fac2*(CONSERVS[STILDEZ] + fac1*Bbar_z);
  return 0;
}

inline int font_fix_general_gamma(CCTK_REAL &u_x, CCTK_REAL &u_y, CCTK_REAL &u_z,CCTK_REAL *CONSERVS,CCTK_REAL *PRIMS,CCTK_REAL *METRIC_PHYS,CCTK_REAL *METRIC_LAP_PSI4,eos_struct &eos) {
  CCTK_REAL gamma=eos.gamma_tab[0];
  CCTK_REAL rhob;
  CCTK_REAL tol = 1.e-15;
  CCTK_REAL gam1 = gamma-1.0;
  CCTK_REAL gam_gamm1_kpoly = gamma/gam1*eos.K_poly;

  CCTK_REAL Bxbar = PRIMS[BX_CENTER]*ONE_OVER_SQRT_4PI;
  CCTK_REAL Bybar = PRIMS[BY_CENTER]*ONE_OVER_SQRT_4PI;
  CCTK_REAL Bzbar = PRIMS[BZ_CENTER]*ONE_OVER_SQRT_4PI;
  CCTK_REAL Bbar_x = METRIC_PHYS[GXX]*Bxbar + METRIC_PHYS[GXY]*Bybar + METRIC_PHYS[GXZ]*Bzbar;
  CCTK_REAL Bbar_y = METRIC_PHYS[GXY]*Bxbar + METRIC_PHYS[GYY]*Bybar + METRIC_PHYS[GYZ]*Bzbar;
  CCTK_REAL Bbar_z = METRIC_PHYS[GXZ]*Bxbar + METRIC_PHYS[GYZ]*Bybar + METRIC_PHYS[GZZ]*Bzbar;
  CCTK_REAL B2bar = Bxbar*Bbar_x + Bybar*Bbar_y + Bzbar*Bbar_z;
  CCTK_REAL Bbar = sqrt(B2bar);

  CCTK_REAL check_B_small = fabs(Bxbar)+fabs(Bybar)+fabs(Bzbar);
  if (check_B_small>0 && check_B_small<1.e-150) {
    // need to compute B2bar specially to prevent floating-point underflow
    CCTK_REAL Bmax = fabs(Bxbar);
    if (Bmax < fabs(Bybar)) Bmax=fabs(Bybar);
    if (Bmax < fabs(Bzbar)) Bmax=fabs(Bzbar);
    CCTK_REAL Bxtmp=Bxbar/Bmax, Bytemp=Bybar/Bmax, Bztemp=Bzbar/Bmax;
    CCTK_REAL B_xtemp=Bbar_x/Bmax, B_ytemp=Bbar_y/Bmax, B_ztemp=Bbar_z/Bmax;
    Bbar = sqrt(Bxtmp*B_xtemp + Bytemp*B_ytemp + Bztemp*B_ztemp)*Bmax;
  }
  CCTK_REAL BbardotS = Bxbar*CONSERVS[STILDEX] + Bybar*CONSERVS[STILDEY] + Bzbar*CONSERVS[STILDEZ];
  CCTK_REAL BbardotS2 = BbardotS*BbardotS;
  CCTK_REAL hatBbardotS = BbardotS/Bbar;
  if (Bbar<1.e-300) hatBbardotS = 0.0;
  CCTK_REAL Psim6 = 1.0/METRIC_LAP_PSI4[PSI6];

  // Limit hatBbardotS
  //CCTK_REAL max_gammav = 100.0;
  //CCTK_REAL rhob_max = CONSERVS[RHOSTAR]*Psim6;
  //CCTK_REAL hmax = 1.0 + gam_gamm1_kpoly*pow(rhob_max,gam1);
  //CCTK_REAL abs_hatBbardotS_max = sqrt(SQR(max_gammav)-1.0)*CONSERVS[RHOSTAR]*hmax;
  //if (fabs(hatBbardotS) > abs_hatBbardotS_max) {
  //   CCTK_REAL fac_reduce = abs_hatBbardotS_max/fabs(hatBbardotS);
  //   CCTK_REAL hatBbardotS_max = hatBbardotS*fac_reduce;
  //   CCTK_REAL Bbar_inv = 1.0/Bbar;
  //   CCTK_REAL hat_Bbar_x = Bbar_x*Bbar_inv;
  //   CCTK_REAL hat_Bbar_y = Bbar_y*Bbar_inv;
  //   CCTK_REAL hat_Bbar_z = Bbar_z*Bbar_inv;
  //   CCTK_REAL sub_fact = hatBbardotS_max - hatBbardotS;
  //   CONSERVS[STILDEX] += sub_fact*hat_Bbar_x;
  //   CONSERVS[STILDEY] += sub_fact*hat_Bbar_y;
  //   CONSERVS[STILDEZ] += sub_fact*hat_Bbar_z;
  //   hatBbardotS = hatBbardotS_max;
  //   BbardotS *= fac_reduce;
  //   BbardotS2 = BbardotS*BbardotS;
  //}

  CCTK_REAL sdots = METRIC_PHYS[GUPXX]*SQR(CONSERVS[STILDEX]) + METRIC_PHYS[GUPYY]*SQR(CONSERVS[STILDEY]) + METRIC_PHYS[GUPZZ]*SQR(CONSERVS[STILDEZ])
    + 2.0*( METRIC_PHYS[GUPXY]*CONSERVS[STILDEX]*CONSERVS[STILDEY] + METRIC_PHYS[GUPXZ]*CONSERVS[STILDEX]*CONSERVS[STILDEZ]
            + METRIC_PHYS[GUPYZ]*CONSERVS[STILDEY]*CONSERVS[STILDEZ]);

  if (sdots<1.e-300) {
    rhob = CONSERVS[RHOSTAR]*Psim6;
    u_x=0.0; u_y=0.0; u_z=0.0;
    return 0;
  }
  /* This test has some problem.
     if (fabs(BbardotS2 - sdots*B2bar) > 1e-8) {
     CCTK_VInfo(CCTK_THORNSTRING,"(Bbar dot S)^2, Bbar^2 * sdotS, %e %e",SQR(BbardotS),sdots*B2bar);
     CCTK_VInfo(CCTK_THORNSTRING,"Cauchy-Schwartz inequality is violated!");
     }
  */
  // Initial guess for W, S_fluid and rhob
  CCTK_REAL W0 = sqrt( SQR(hatBbardotS) + SQR(CONSERVS[RHOSTAR]) ) * Psim6;
  CCTK_REAL Sf20 = (SQR(W0)*sdots + BbardotS2*(B2bar + 2.0*W0))/SQR(W0+B2bar);
  CCTK_REAL rhob0 = CONSERVS[RHOSTAR]*Psim6/sqrt(1.0+Sf20/SQR(CONSERVS[RHOSTAR]));
  CCTK_REAL W=W0,Sf2=Sf20,rhob1=rhob0;

  //****************************************************************
  //                          FONT FIX
  // Impose Font fix when HARM primitives solver fails to find
  //   acceptable set of primitives.
  //****************************************************************

  bool fontcheck=true;

  int itcount = 0, maxits=500;
  while(fontcheck && itcount < maxits) {
    itcount++;
    W0 = W;
    Sf20 = Sf2;
    rhob0 = rhob1;

    // first, find rhob for the given S_fluid^2
    CCTK_REAL h = 1.0 + gam_gamm1_kpoly*pow(rhob0,gam1); 
    rhob1 = CONSERVS[RHOSTAR]*Psim6/sqrt(1.0+Sf20/SQR(CONSERVS[RHOSTAR]*h));
    while( fabs(rhob1-rhob0) > rhob1*tol) {
      rhob0 = rhob1;
      h = 1.0 + gam_gamm1_kpoly*pow(rhob0,gam1);
      rhob1 = CONSERVS[RHOSTAR]*Psim6/sqrt(1.0+Sf20/SQR(CONSERVS[RHOSTAR]*h));
    }

    h = 1.0 + gam_gamm1_kpoly*pow(rhob1,gam1);
    W = sqrt( Sf20 + SQR(CONSERVS[RHOSTAR]*h))*Psim6;
    Sf2 = (SQR(W)*sdots + BbardotS2*(B2bar + 2.0*W))/SQR(W+B2bar);
    if ( fabs(W-W0) < W*tol && fabs(Sf20-Sf2) < Sf2*tol) fontcheck=false;
  }

  if (itcount>=maxits) {
    // Increase tol and try again
    fontcheck=true;
    tol *=100.0;
    itcount = 0;
    while(fontcheck && itcount < maxits) {
      itcount++;
      W0 = W;
      Sf20 = Sf2;
      rhob0 = rhob1;

      // first, find rhob for the given S_fluid^2
      CCTK_REAL h = 1.0 + gam_gamm1_kpoly*pow(rhob0,gam1);
      rhob1 = CONSERVS[RHOSTAR]*Psim6/sqrt(1.0+Sf20/SQR(CONSERVS[RHOSTAR]*h));
      while( fabs(rhob1-rhob0) > rhob1*tol) {
        rhob0 = rhob1;
        h = 1.0 + gam_gamm1_kpoly*pow(rhob0,gam1);
        rhob1 = CONSERVS[RHOSTAR]*Psim6/sqrt(1.0+Sf20/SQR(CONSERVS[RHOSTAR]*h));
      }

      h = 1.0 + gam_gamm1_kpoly*pow(rhob1,gam1);
      W = sqrt( Sf20 + SQR(CONSERVS[RHOSTAR]*h))*Psim6;
      Sf2 = (SQR(W)*sdots + BbardotS2*(B2bar + 2.0*W))/SQR(W+B2bar);
      if ( fabs(W-W0) < W*tol && fabs(Sf20-Sf2) < Sf2*tol) fontcheck=false;
    }
  }
  //************************************************************************************************************** 

  if(fontcheck==true) {
    return 1;
  }

  // Font fix works, now compute u_i
  rhob = rhob1;
  CCTK_REAL h = 1.0+gam_gamm1_kpoly*pow(rhob,gam1);
  CCTK_REAL gammav = CONSERVS[RHOSTAR]*Psim6/rhob;
  CCTK_REAL rhosh = CONSERVS[RHOSTAR]*h;
  CCTK_REAL fac1 = METRIC_LAP_PSI4[PSI6]*BbardotS/(gammav*rhosh);
  CCTK_REAL fac2 = 1.0/(rhosh + METRIC_LAP_PSI4[PSI6]*B2bar/gammav);
  u_x = fac2*(CONSERVS[STILDEX] + fac1*Bbar_x);
  u_y = fac2*(CONSERVS[STILDEY] + fac1*Bbar_y);
  u_z = fac2*(CONSERVS[STILDEZ] + fac1*Bbar_z);

  return 0;
}
