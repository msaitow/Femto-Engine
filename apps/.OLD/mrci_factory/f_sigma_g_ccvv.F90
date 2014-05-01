#include "../f_ct.fh"


!      ______                  __           
!     / ____/___   ____ ___   / /_ ____     
!    / /_   / _ \ / __ `__ \ / __// __ \ 
!   / __/  /  __// / / / / // /_ / /_/ /    
!  /_/     \___//_/ /_/ /_/ \__/ \____/  



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ccvv_no0_x0(sv1, iv1, T2, V2, S0, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_sigma_g_ccvv_no0_x0(sv1, iv1, av2_i, h2_i, S0, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_g_ccvv_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_ccvv_no0_x0(s_v1, i_v1, T2_, V2_, S0_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: S0_
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_c1, i_c1, s_v2, i_v2
! S0()  <-- 
! (    4.00000000) T2(c2,c1,v2,v1) V2(v1,c1,c2,v2) 
do s_c2 = 0, nir-1
do s_c1 = 0, nir-1
do s_v2 = 0, nir-1
if( &
IEOR(s_c2,s_c1) == IEOR(s_v2,s_v1) .and. &
IEOR(s_v1,s_c1) == IEOR(s_c2,s_v2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_v2 = psym(I_BEGIN, I_V, s_v2), psym(I_END, I_V, s_v2)
S0_ = S0_ &
  + 4.00000000d+00 & 
  * T2_(s_v2, s_c1, s_c2)%array(i_v2, i_c1, i_c2) & 
  * V2_(s_v2, s_c2, s_c1)%array(i_v2, i_c2, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_g_ccvv_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ccvv_no0_x1(sv2, iv2, T2, V2, S0, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv2, iv2
real(kind=8), intent(inout) :: T2(*), V2(*), S0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv2, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv2, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_sigma_g_ccvv_no0_x1(sv2, iv2, av2_i, h2_i, S0, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_g_ccvv_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_ccvv_no0_x1(s_v2, i_v2, T2_, V2_, S0_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v2, s_v2
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: S0_
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_c2, i_c2, s_v1, i_v1
! S0()  <-- 
! (   -2.00000000) T2(c1,c2,v1,v2) V2(v2,c1,c2,v1) 
do s_c1 = 0, nir-1
do s_c2 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c1,s_c2) == IEOR(s_v1,s_v2) .and. &
IEOR(s_v2,s_c1) == IEOR(s_c2,s_v1)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S0_ = S0_ &
  - 2.00000000d+00 & 
  * T2_(s_v1, s_c2, s_c1)%array(i_v1, i_c2, i_c1) & 
  * V2_(s_v1, s_c2, s_c1)%array(i_v1, i_c2, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_g_ccvv_no0_x1

