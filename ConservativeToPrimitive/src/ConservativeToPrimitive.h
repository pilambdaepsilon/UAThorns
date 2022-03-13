#ifndef CONSERVATIVETOPRIMITIVE_H_
#define CONSERVATIVETOPRIMITIVE_H_

#define C2P_RHO 0
#define C2P_v1_cov 1
#define C2P_v2_cov 2
#define C2P_v3_cov 3
#define C2P_EPS 4

#define C2P_D 0
#define C2P_S1_cov 1
#define C2P_S2_cov 2
#define C2P_S3_cov 3
#define C2P_TAU 4

#define C2P_B1_con 5
#define C2P_B2_con 6
#define C2P_B3_con 7

#define C2P_YE 8
#define C2P_TEMP 9
#define C2P_PRESS 10
#define C2P_WLORENTZ 11
#define C2P_ENT 12
#define C2P_A_BAR 13
#define C2P_MUHAT 14

static const int C2P_PHI=0,C2P_PSI=1,C2P_GXX=2,C2P_GXY=3,C2P_GXZ=4,C2P_GYY=5,C2P_GYZ=6,C2P_GZZ=7,
  C2P_LAPM1=8,C2P_SHIFTX=9,C2P_SHIFTY=10,C2P_SHIFTZ=11,C2P_GUPXX=12,C2P_GUPYY=13,C2P_GUPZZ=14;
static const int C2P_GUPXY=15,C2P_GUPXZ=16,C2P_GUPYZ=17;
static const int C2P_LAPSE=0,C2P_PSI2=1,C2P_PSI4=2,C2P_PSI6=3,C2P_PSIM4=4,C2P_LAPSEINV=5;
#endif
