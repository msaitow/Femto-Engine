#include "../f_ct.fh"


!  `MMMMMMM                                         
!   MM    \                         /              
!   MM       ____  ___  __    __   /M      _____    
!   MM   ,  6MMMMb `MM 6MMb  6MMb /MMMMM  6MMMMMb   
!   MMMMMM 6M'  `Mb MM69 `MM69 `Mb MM    6M'   `Mb  
!   MM   ` MM    MM MM'   MM'   MM MM    MM     MM  
!   MM     MMMMMMMM MM    MM    MM MM    MM     MM  
!   MM     MM       MM    MM    MM MM    MM     MM  
!   MM     YM    d9 MM    MM    MM YM.  ,YM.   ,M9  
!  _MM_     YMMMM9 _MM_  _MM_  _MM_ YMMM9 YMMMMM9   



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ooov_no0_x0(sc, ic, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc

call set_symblock_xccaaa(sleft, x, nir, nsym, psym) ! -> xccaaa (allocate) 
call g_sigma_ccvv_ooov_no0_x0(sc, ic, h2_i, xccaaa, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(xccaaa)

end subroutine g_if_sigma_ccvv_ooov_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ooov_no0_x0(s_c, i_c, V2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3, s_o2, i_o2, s_o4, i_o4, s_w, i_w
integer :: s_y, i_y
! X(w,y,o1,o2,o4,c)  <-- 
! (    1.00000000)  D2(o1,o3,o2,o4) V2(c,w,y,o3) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(IEOR(s_w,s_y),s_o1) == IEOR(IEOR(s_o2,s_o4),s_c) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_c,s_w) == IEOR(s_y,s_o3)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
X_(s_o4, s_o2, s_o1, s_y, s_w)%array(i_o4, i_o2, i_o1, i_y, i_w) =  &
    X_(s_o4, s_o2, s_o1, s_y, s_w)%array(i_o4, i_o2, i_o1, i_y, i_w) &
  + 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * V2_(s_o3, s_y, s_w)%array(i_o3, i_y, i_w)
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

end subroutine g_sigma_ccvv_ooov_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ooov_no1_x0(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc

call set_symblock_xccaaa(sleft, x, nir, nsym, psym) ! -> xccaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ooov_no1_x0(sa, ia, sc, ic, av2_i, xccaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xccaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ooov_no1_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ooov_no1_x0(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_o4, i_o4, s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    1.00000000) T2(o2,o1,o4,a) X(w,y,o1,o2,o4,c) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_o4,s_a) .and. &
IEOR(IEOR(s_w,s_y),s_o1) == IEOR(IEOR(s_o2,s_o4),s_c)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * T2_(s_o4, s_o1, s_o2)%array(i_o4, i_o1, i_o2) & 
  * X_(s_o4, s_o2, s_o1, s_y, s_w)%array(i_o4, i_o2, i_o1, i_y, i_w)
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

end subroutine g_sigma_ccvv_ooov_no1_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ooov_no0_x1(sc, ic, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc

call set_symblock_xccaaa(sleft, x, nir, nsym, psym) ! -> xccaaa (allocate) 
call g_sigma_ccvv_ooov_no0_x1(sc, ic, h2_i, xccaaa, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(xccaaa)

end subroutine g_if_sigma_ccvv_ooov_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ooov_no0_x1(s_c, i_c, V2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3, s_o2, i_o2, s_o4, i_o4, s_y, i_y
integer :: s_w, i_w
! X(y,w,o1,o2,o4,c)  <-- 
! (    1.00000000)  D2(o1,o3,o2,o4) V2(c,y,w,o3) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(IEOR(s_y,s_w),s_o1) == IEOR(IEOR(s_o2,s_o4),s_c) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_c,s_y) == IEOR(s_w,s_o3)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
X_(s_o4, s_o2, s_o1, s_w, s_y)%array(i_o4, i_o2, i_o1, i_w, i_y) =  &
    X_(s_o4, s_o2, s_o1, s_w, s_y)%array(i_o4, i_o2, i_o1, i_w, i_y) &
  + 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * V2_(s_o3, s_w, s_y)%array(i_o3, i_w, i_y)
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

end subroutine g_sigma_ccvv_ooov_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ooov_no1_x1(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc

call set_symblock_xccaaa(sleft, x, nir, nsym, psym) ! -> xccaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ooov_no1_x1(sa, ia, sc, ic, av2_i, xccaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xccaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ooov_no1_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ooov_no1_x1(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_o4, i_o4, s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(o2,o1,o4,a) X(y,w,o1,o2,o4,c) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_o4,s_a) .and. &
IEOR(IEOR(s_y,s_w),s_o1) == IEOR(IEOR(s_o2,s_o4),s_c)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_o4, s_o1, s_o2)%array(i_o4, i_o1, i_o2) & 
  * X_(s_o4, s_o2, s_o1, s_w, s_y)%array(i_o4, i_o2, i_o1, i_w, i_y)
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

end subroutine g_sigma_ccvv_ooov_no1_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ooov_no0_x2(sa, ia, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xccaaa(sleft, x, nir, nsym, psym) ! -> xccaaa (allocate) 
call g_sigma_ccvv_ooov_no0_x2(sa, ia, h2_i, xccaaa, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(xccaaa)

end subroutine g_if_sigma_ccvv_ooov_no0_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ooov_no0_x2(s_a, i_a, V2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3, s_o2, i_o2, s_o4, i_o4, s_y, i_y
integer :: s_w, i_w
! X(y,w,o1,o2,o4,a)  <-- 
! (    1.00000000)  D2(o1,o3,o2,o4) V2(a,y,w,o3) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(IEOR(s_y,s_w),s_o1) == IEOR(IEOR(s_o2,s_o4),s_a) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_a,s_y) == IEOR(s_w,s_o3)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
X_(s_o4, s_o2, s_o1, s_w, s_y)%array(i_o4, i_o2, i_o1, i_w, i_y) =  &
    X_(s_o4, s_o2, s_o1, s_w, s_y)%array(i_o4, i_o2, i_o1, i_w, i_y) &
  + 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * V2_(s_o3, s_w, s_y)%array(i_o3, i_w, i_y)
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

end subroutine g_sigma_ccvv_ooov_no0_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ooov_no1_x2(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_xccaaa(sleft, x, nir, nsym, psym) ! -> xccaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ooov_no1_x2(sa, ia, sc, ic, av2_i, xccaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xccaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ooov_no1_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ooov_no1_x2(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_o4, i_o4, s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    1.00000000) T2(o2,o1,o4,c) X(y,w,o1,o2,o4,a) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_o4,s_c) .and. &
IEOR(IEOR(s_y,s_w),s_o1) == IEOR(IEOR(s_o2,s_o4),s_a)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * T2_(s_o4, s_o1, s_o2)%array(i_o4, i_o1, i_o2) & 
  * X_(s_o4, s_o2, s_o1, s_w, s_y)%array(i_o4, i_o2, i_o1, i_w, i_y)
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

end subroutine g_sigma_ccvv_ooov_no1_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ooov_no0_x3(sa, ia, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xccaaa(sleft, x, nir, nsym, psym) ! -> xccaaa (allocate) 
call g_sigma_ccvv_ooov_no0_x3(sa, ia, h2_i, xccaaa, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(xccaaa)

end subroutine g_if_sigma_ccvv_ooov_no0_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ooov_no0_x3(s_a, i_a, V2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3, s_o2, i_o2, s_o4, i_o4, s_w, i_w
integer :: s_y, i_y
! X(w,y,o1,o3,o4,a)  <-- 
! (    1.00000000)  D2(o1,o3,o2,o4) V2(a,w,y,o2) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(IEOR(s_w,s_y),s_o1) == IEOR(IEOR(s_o3,s_o4),s_a) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_a,s_w) == IEOR(s_y,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
X_(s_o4, s_o3, s_o1, s_y, s_w)%array(i_o4, i_o3, i_o1, i_y, i_w) =  &
    X_(s_o4, s_o3, s_o1, s_y, s_w)%array(i_o4, i_o3, i_o1, i_y, i_w) &
  + 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * V2_(s_o2, s_y, s_w)%array(i_o2, i_y, i_w)
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

end subroutine g_sigma_ccvv_ooov_no0_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ooov_no1_x3(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_xccaaa(sleft, x, nir, nsym, psym) ! -> xccaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ooov_no1_x3(sa, ia, sc, ic, av2_i, xccaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xccaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ooov_no1_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ooov_no1_x3(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o4, i_o4, s_o1, i_o1, s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(o3,o4,o1,c) X(w,y,o1,o3,o4,a) 
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o3,s_o4) == IEOR(s_o1,s_c) .and. &
IEOR(IEOR(s_w,s_y),s_o1) == IEOR(IEOR(s_o3,s_o4),s_a)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_o1, s_o4, s_o3)%array(i_o1, i_o4, i_o3) & 
  * X_(s_o4, s_o3, s_o1, s_y, s_w)%array(i_o4, i_o3, i_o1, i_y, i_w)
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

end subroutine g_sigma_ccvv_ooov_no1_x3

