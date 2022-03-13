#ifndef LOOP_DEFINES_RECONSTRUCTION_H_
#define LOOP_DEFINES_RECONSTRUCTION_H_

#define LOOP_DEFINE(gz_shift_lo,gz_shift_hi,  ext,flux_dirn,  ijkgz_lo_hi,gz_lo,gz_hi) \
  for(int rr=1;rr<=3;rr++) {                                            \
    ijkgz_lo_hi[rr][0]=          gz_lo[rr];                             \
    ijkgz_lo_hi[rr][1]=ext[rr-1]-gz_hi[rr];                             \
  }                                                                     \
  ijkgz_lo_hi[flux_dirn][0] += gz_shift_lo;                             \
  ijkgz_lo_hi[flux_dirn][1] -= gz_shift_hi;                             \
  /* The following line is valid C99 */                                 \
  _Pragma("omp parallel for private(U,dU,slope_lim_dU,Ur,Ul)")          \
  for(int k=ijkgz_lo_hi[3][0];k<ijkgz_lo_hi[3][1];k++)                  \
    for(int j=ijkgz_lo_hi[2][0];j<ijkgz_lo_hi[2][1];j++)                \
      for(int i=ijkgz_lo_hi[1][0];i<ijkgz_lo_hi[1][1];i++)

// This define only sets indices. 
// FIXME: benchmark with and without the if() statement.
// FIXME: try without index_arr being defined in all directions.
#define SET_INDEX_ARRAYS(IMIN,IMAX,flux_dirn)                           \
  int max_shift=(MAXNUMINDICES/2);                                      \
  /* DEBUGGING ONLY:  if(IMIN<-max_shift || IMAX>max_shift) CCTK_VError(VERR_DEF_PARAMS,"FIX MAXNUMINDICES!"); */ \
  int index_arr[4][MAXNUMINDICES];                                      \
  for(int idx=IMIN;idx<=IMAX;idx++) {                                   \
    index_arr[flux_dirn][idx+max_shift]=                                \
      CCTK_GFINDEX3D(cctkGH,                                            \
                     i+idx*kronecker_delta[flux_dirn][0],               \
                     j+idx*kronecker_delta[flux_dirn][1],               \
                     k+idx*kronecker_delta[flux_dirn][2]);              \
  }

#define SET_INDEX_ARRAYS_3DBLOCK(IJKLOHI)                               \
  int max_shift=(MAXNUMINDICES/2);                                      \
  int index_arr_3DB[MAXNUMINDICES][MAXNUMINDICES][MAXNUMINDICES];       \
  for(int idx_k=IJKLOHI[4];idx_k<=IJKLOHI[5];idx_k++) for(int idx_j=IJKLOHI[2];idx_j<=IJKLOHI[3];idx_j++) for(int idx_i=IJKLOHI[0];idx_i<=IJKLOHI[1];idx_i++) { \
        index_arr_3DB[idx_k+max_shift][idx_j+max_shift][idx_i+max_shift]=CCTK_GFINDEX3D(cctkGH,i+idx_i,j+idx_j,k+idx_k); \
      }

#endif /* LOOP_DEFINES_RECONSTRUCTION_H_ */
