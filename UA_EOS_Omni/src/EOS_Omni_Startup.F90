#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine EOS_Omni_Startup(CCTK_ARGUMENTS)

  use EOS_Omni_Module
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  integer :: p
  
  if(poly_gamma_initial .gt. 0d0) then
    poly_gamma_ini = poly_gamma_initial
  else
    poly_gamma_ini = poly_gamma
  end if

  poly_k_cgs = poly_k * rho_gf**poly_gamma_ini / press_gf

  gl_k_cgs   = gl_k * rho_gf**poly_gamma_ini / press_gf

! --------------------------------------------
! Use Legacy parameter in case n_pieces=2 
!       and hybrid_gamma(1) has not been set
! --------------------------------------------
    
  if(n_pieces .gt. 0) then
     allocate(hybrid_eps(n_pieces))
     allocate(hybrid_k(n_pieces))
     hybrid_k(1) = hybrid_k0
     hybrid_eps(1) = 0
     if (n_pieces .gt. 1) then
        do p = 1,n_pieces-1
           hybrid_k(p+1) = hybrid_k(p) * hybrid_rho(p)**(hybrid_gamma(p)-hybrid_gamma(p+1))
           hybrid_eps(p+1) = hybrid_eps(p) + &
                hybrid_k(p)*hybrid_rho(p)**(hybrid_gamma(p)-1.0d0)/(hybrid_gamma(p)-1.0d0) - &
                hybrid_k(p+1) * hybrid_rho(p)**(hybrid_gamma(p+1)-1.0d0)/(hybrid_gamma(p+1)-1.0d0) 
        end do
     end if
  end if
  if( (n_pieces .eq. 2) .and. (hybrid_gamma(1) .eq. 0d0 ) )then
     hybrid_k(1)     = hybrid_k1
     if (poly_gamma_initial .gt. 0d0) then
       hybrid_k(1)   = hybrid_k(1) * rho_gf**(hybrid_gamma(1) - poly_gamma_initial)
     end if
     hybrid_rho(1)   = hybrid_rho_nuc
     hybrid_gamma(1) = hybrid_gamma1
     hybrid_gamma(2) = hybrid_gamma2
  end if


end subroutine EOS_Omni_Startup

subroutine EOS_Omni_Get_Energy_Shift(CCTK_ARGUMENTS)

  use EOS_Omni_Module
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

!  call nuc_eos_c_get_energy_shift(energy_shift,eos_tempmin,eos_tempmax,&
!       eos_yemin,eos_yemax)
  call nuc_eos_c_get_energy_shift(energy_shift,eos_tempmin,eos_tempmax,&
       eos_yemin,eos_yemax,eos_rhomin,eos_rhomax,eos_epsmin,eos_epsmax)

end subroutine EOS_Omni_Get_Energy_Shift
