//--------------------------------------------------------------------------
// Register with the time stepper 
// (MoL thorn, found in arrangements/CactusBase/MoL)
// To understand this, read documentation in arrangements/CactusBase/MoL/doc
//--------------------------------------------------------------------------

#include "cctk.h"
#include <cstdio>
#include <cmath>
#include <cstddef>
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

extern "C" void IllinoisGRMHD_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0, group, rhs;

  // Register evolution & RHS gridfunction groups

  /* Ax and Ax_rhs */
  group = CCTK_GroupIndex("IllinoisGRMHD::em_Ax");
  rhs = CCTK_GroupIndex("IllinoisGRMHD::em_Ax_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Ay and Ay_rhs */
  group = CCTK_GroupIndex("IllinoisGRMHD::em_Ay");
  rhs = CCTK_GroupIndex("IllinoisGRMHD::em_Ay_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Az and Az_rhs */
  group = CCTK_GroupIndex("IllinoisGRMHD::em_Az");
  rhs = CCTK_GroupIndex("IllinoisGRMHD::em_Az_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* psi6phi and psi6phi_rhs */
  group = CCTK_GroupIndex("IllinoisGRMHD::em_psi6phi");
  rhs = CCTK_GroupIndex("IllinoisGRMHD::em_psi6phi_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* ALL OTHER EVOLVED VARIABLES (rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z) */
  group = CCTK_GroupIndex("IllinoisGRMHD::grmhd_conservatives");
  rhs = CCTK_GroupIndex("IllinoisGRMHD::grmhd_conservatives_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");

}
