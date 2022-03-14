//
// This subroutine calcualtes the eigenvalues of a real, symmetric 3x3
// matrix M={{M11,M12,M13},{M12,M22,M23},{M13,M23,M33}} based on the
// algorithm described in
// http://en.wikipedia.org/wiki/Eigenvalue_algorithm#Eigenvalues_of_3.C3.973_matrices
// which simply solve the cubic equation Det( M - lamnda I)=0 analytically.
// The eigenvalues are stored in lam1, lam2 and lam3.
//
void eigenvalues_3by3_real_sym_matrix(CCTK_REAL & lam1, CCTK_REAL & lam2, CCTK_REAL & lam3,
                                      CCTK_REAL M11, CCTK_REAL M12, CCTK_REAL M13, CCTK_REAL M22, CCTK_REAL M23, CCTK_REAL M33)
{
  CCTK_REAL m = (M11 + M22 + M33)/3.0;
  CCTK_REAL K11 = M11 - m, K12 = M12, K13 = M13, K22 = M22-m, K23 = M23, K33=M33-m;
  CCTK_REAL q = 0.5* (K11*K22*K33 + K12*K23*K13 + K13*K12*K23 - K13*K22*K13
                      - K12*K12*K33 - K11*K23*K23);
  CCTK_REAL p = ( SQR(K11) + SQR(K22) + SQR(K33) + 2.0*(SQR(K12) + SQR(K13) + SQR(K23) ) )/6.0;
  CCTK_REAL phi;
  CCTK_REAL p32 = sqrt(p*p*p);
  if (fabs(q) >= fabs(p32) ) {
    phi = 0.0;
  } else {
    phi = acos(q/p32)/3.0;
  }
  if (phi<0.0) phi += M_PI/3.0;
  CCTK_REAL sqrtp = sqrt(p);
  CCTK_REAL sqrtp_cosphi = sqrtp*cos(phi);
  CCTK_REAL sqrtp_sqrt3_sinphi = sqrtp*sqrt(3.0)*sin(phi);
  lam1 = m + 2.0*sqrtp_cosphi;
  lam2 = m - sqrtp_cosphi - sqrtp_sqrt3_sinphi;
  lam3 = m - sqrtp_cosphi + sqrtp_sqrt3_sinphi;
}
