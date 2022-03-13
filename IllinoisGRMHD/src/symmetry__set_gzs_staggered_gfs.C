#include "cctk.h"
#include "cctk_Parameters.h"
#include <cstdio>
#include <cstdlib>
#include "IllinoisGRMHD_headers.h"

void IllinoisGRMHD_set_symmetry_gzs_staggered(const cGH *cctkGH, const int *cctk_lsh,CCTK_REAL *X,CCTK_REAL *Y,CCTK_REAL *Z, CCTK_REAL *gridfunc,
                                              CCTK_REAL *gridfunc_syms,int stagger_x,int stagger_y,int stagger_z) {

  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(Symmetry, "equatorial")) 
    CCTK_VError(VERR_DEF_PARAMS,"Warning: Symmetry==equatorial not supported! USE AT YOUR OWN RISK. You will need to comment this error message out.");

  // No symmetries -> return.
  if(CCTK_EQUALS(Symmetry, "none")) return;

  CCTK_REAL dz = Z[CCTK_GFINDEX3D(cctkGH,0,0,1)] - Z[CCTK_GFINDEX3D(cctkGH,0,0,0)];

  CCTK_REAL z_offset = dz*0.5*stagger_z;

  int num_gzs=0;
  //FIXME: Might want to use cctk_nghostzones instead...
  while( (Z[CCTK_GFINDEX3D(cctkGH,0,0,num_gzs)]+z_offset) < -dz*0.1 && num_gzs<cctk_lsh[2]) num_gzs++;
  if(num_gzs*2>=cctk_lsh[2]) CCTK_VError(VERR_DEF_PARAMS,"ERROR in symmetry__set_gzs_staggered_gfs.C");

#pragma omp parallel for
  for(int k=0;k<num_gzs;k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        int index_inside__sym_gz = CCTK_GFINDEX3D(cctkGH,i,j,k);

        /* This loop sets symmetry ghostzones, regardless of how the gridfunction is staggered.
         *
         * STAGGERED PATTERN:
         * if num_gzs==1 && stagger_z==1: 
         * z[] = {-dz/2,dz/2,3dz/2, etc} -> gridfunc[index 0] = gridfunc_syms[2]*gridfunc[index 1]
         *
         * if num_gzs==2 && stagger_z==1: 
         * z[] = {-3dz/2,-dz/2,dz/2,3dz/2 etc} 
         * -> gridfunc[index 0] = gridfunc_syms[2]*gridfunc[index 3]
         * -> gridfunc[index 1] = gridfunc_syms[2]*gridfunc[index 2]
         * .
         * .
         * .
         * -> gridfunc[i] = gridfunc_syms[2]*gridfunc[(num_gz*2-1)-i]
         *
         * UNSTAGGERED PATTERN:
         * if num_gzs==1 && stagger_z==0:
         * z[] = {-dz,0,dz, etc} -> gridfunc[index 0] = gridfunc_syms[2]*gridfunc[index 2]
         *
         * if num_gzs==2 && stagger_z==0:
         * z[] = {-2dz,-dz,0,dz,2dz, etc} -> gridfunc[index 0] = gridfunc_syms[2]*gridfunc[index 4]
         * z[] = {-2dz,-dz,0,dz,2dz, etc} -> gridfunc[index 1] = gridfunc_syms[2]*gridfunc[index 3]
         * .
         * .
         * .
         * -> gridfunc[i] = gridfunc_syms[2]*gridfunc[(num_gz*2)-i]
         *
         * OVERALL PATTERN: gridfunc[i] = gridfunc_syms[2]*gridfunc[(num_gz*2-stagger_z)-i] */

        int matching_index_outside_sym_gz = CCTK_GFINDEX3D(cctkGH,i,j,(num_gzs*2-stagger_z)-k);

        gridfunc[index_inside__sym_gz] = gridfunc_syms[2]*gridfunc[matching_index_outside_sym_gz];
      }
}
