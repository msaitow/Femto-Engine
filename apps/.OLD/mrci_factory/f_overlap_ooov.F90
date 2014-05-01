#include "../f_ct.fh"


!  __/\\\\\\\\\\\\\\\____________________________________________________________________                                   
!   _\/\\\///////////_____________________________________________________________________                                             
!    _\/\\\_______________________________________________________/\\\_____________________                                         
!     _\/\\\\\\\\\\\__________/\\\\\\\\______/\\\\\__/\\\\\_____/\\\\\\\\\\\______/\\\\\____ 
!      _\/\\\///////_________/\\\/////\\\___/\\\///\\\\\///\\\__\////\\\////_____/\\\///\\\__               
!       _\/\\\_______________/\\\\\\\\\\\___\/\\\_\//\\\__\/\\\_____\/\\\________/\\\__\//\\\_       
!        _\/\\\______________\//\\///////____\/\\\__\/\\\__\/\\\_____\/\\\_/\\___\//\\\__/\\\__            
!         _\/\\\_______________\//\\\\\\\\\\__\/\\\__\/\\\__\/\\\_____\//\\\\\_____\///\\\\\/___    
!          _\///_________________\//////////___\///___\///___\///_______\/////________\/////_____                                   



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_overlap_ooov_no0_x0(sa, ia, T2, O2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), O2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, O2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_overlap_ooov_no0_x0(sa, ia, av2_i, av2_i2, d3, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_overlap_ooov_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_overlap_ooov_no0_x0(s_a, i_a, T2_, O2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: O2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_m, i_m, s_o2, i_o2
integer :: s_o3, i_o3
! O2(i,k,m,a)  <-- 
! (    1.00000000) D3(i,m,k,o2,o1,o3) T2(o3,o2,o1,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o2,s_o1),s_o3) .and. &
IEOR(s_o3,s_o2) == IEOR(s_o1,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
O2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    O2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o3, s_o1, s_o2, s_k, s_m, s_i)%array(i_o3, i_o1, i_o2, i_k, i_m, i_i) & 
  * T2_(s_o1, s_o2, s_o3)%array(i_o1, i_o2, i_o3)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_overlap_ooov_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_overlap_ooov_no0_x1(sa, ia, T2, O2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), O2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, O2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_overlap_ooov_no0_x1(sa, ia, av2_i, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_overlap_ooov_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_overlap_ooov_no0_x1(s_a, i_a, T2_, O2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: O2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2, s_m, i_m
! O2(i,k,m,a)  <-- 
! (    1.00000000) D2(i,o1,k,o2) T2(o1,o2,m,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_o2) .and. &
IEOR(s_o1,s_o2) == IEOR(s_m,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
O2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    O2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i) & 
  * T2_(s_m, s_o2, s_o1)%array(i_m, i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_overlap_ooov_no0_x1

