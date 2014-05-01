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
subroutine g_if_overlap_ccvv_no0_x0(sc, ic, T2, O2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), O2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, O2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_overlap_ccvv_no0_x0(sc, ic, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_overlap_ccvv_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_overlap_ccvv_no0_x0(s_c, i_c, T2_, O2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: O2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! O2(w,y,a,c)  <-- 
! (    4.00000000) T2(w,y,a,c) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_a,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
O2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    O2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_a, s_y, s_w)%array(i_a, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_overlap_ccvv_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_overlap_ccvv_no0_x1(sc, ic, T2, O2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), O2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, O2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_overlap_ccvv_no0_x1(sc, ic, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_overlap_ccvv_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_overlap_ccvv_no0_x1(s_c, i_c, T2_, O2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: O2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! O2(w,y,a,c)  <-- 
! (   -2.00000000) T2(y,w,a,c) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_a,s_c)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
O2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    O2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_a, s_w, s_y)%array(i_a, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_overlap_ccvv_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_overlap_ccvv_no0_x2(sa, ia, sc, ic, T2, O2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), O2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, O2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_overlap_ccvv_no0_x2(sa, ia, sc, ic, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_overlap_ccvv_no0_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_overlap_ccvv_no0_x2(s_a, i_a, s_c, i_c, T2_, O2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: O2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! O2(w,y,a,c)  <-- 
! (    4.00000000) T2(y,w,c,a) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_c,s_a)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
O2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    O2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_c, s_w, s_y)%array(i_c, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_overlap_ccvv_no0_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_overlap_ccvv_no0_x3(sa, ia, sc, ic, T2, O2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), O2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, O2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_overlap_ccvv_no0_x3(sa, ia, sc, ic, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_overlap_ccvv_no0_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_overlap_ccvv_no0_x3(s_a, i_a, s_c, i_c, T2_, O2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: O2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! O2(w,y,a,c)  <-- 
! (   -2.00000000) T2(w,y,c,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_c,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
O2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    O2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_c, s_y, s_w)%array(i_c, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_overlap_ccvv_no0_x3

