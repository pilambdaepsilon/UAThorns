#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine driver_reset_gamma(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  integer:: i

  write(*,*) "driver_reset_gamma,dx=",CCTK_DELTA_SPACE(1)

  if (RESET_GAMMA==1) then
     gamma_th = gamma_th_new
     eps_tab(1)=K_poly/(gamma_th-1.0d0)
     do i=1,neos
        gamma_tab(i)=gamma_th_new
     enddo
  end if

  write(*,*) "new Gamma_th = ",gamma_th


end subroutine driver_reset_gamma
