#include "../f_ct.fh"


!      ______                  __           
!     / ____/___   ____ ___   / /_ ____     
!    / /_   / _ \ / __ `__ \ / __// __ \ 
!   / __/  /  __// / / / / // /_ / /_/ /    
!  /_/     \___//_/ /_/ /_/ \__/ \____/  



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x0(sa, ia, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x0(sa, ia, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x0(s_a, i_a, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o3, i_o3, s_o1, i_o1, s_o4, i_o4
! X(o2,o3,o4,a)  <-- 
! (    1.00000000)  T2(o2,o3,o1,a) h(o4,o1) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o2,s_o3) == IEOR(s_o4,s_a) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_o1,s_a) .and. &
IEOR(s_o4,s_o1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o3, s_o2)%array(i_o4, i_o3, i_o2) =  &
    X_(s_o4, s_o3, s_o2)%array(i_o4, i_o3, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2) & 
  * h_(s_o1, s_o4)%array(i_o1, i_o4)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x0(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x0(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x0(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o3, i_o3, s_o4, i_o4
integer :: s_o2, i_o2
! S2(i,k,m,a)  <-- 
! (    1.00000000) D3(i,m,k,o3,o4,o2) X(o2,o3,o4,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o3,s_o4),s_o2) .and. &
IEOR(s_o2,s_o3) == IEOR(s_o4,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o4, s_o3, s_k, s_m, s_i)%array(i_o2, i_o4, i_o3, i_k, i_m, i_i) & 
  * X_(s_o4, s_o3, s_o2)%array(i_o4, i_o3, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x1(sa, ia, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x1(sa, ia, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x1(s_a, i_a, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_o1, i_o1, s_m, i_m
! X(o3,o2,m,a)  <-- 
! (    1.00000000)  T2(o3,o2,o1,a) h(m,o1) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o3,s_o2) == IEOR(s_m,s_a) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_o1,s_a) .and. &
IEOR(s_m,s_o1) == 0) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o2, s_o3)%array(i_m, i_o2, i_o3) =  &
    X_(s_m, s_o2, s_o3)%array(i_m, i_o2, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o2, s_o3)%array(i_o1, i_o2, i_o3) & 
  * h_(s_o1, s_m)%array(i_o1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x1(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x1(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x1(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_m, i_m
! S2(i,k,m,a)  <-- 
! (    1.00000000) D2(i,o3,k,o2) X(o3,o2,m,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o3,s_o2) == IEOR(s_m,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * X_(s_m, s_o2, s_o3)%array(i_m, i_o2, i_o3)
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

end subroutine g_sigma_ooov_ooov_no1_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x2(sa, ia, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x2(sa, ia, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x2(s_a, i_a, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3, s_o4, i_o4, s_o2, i_o2
! X(o3,o4,o2,a)  <-- 
! (    1.00000000)  T2(o1,o3,o4,a) h(o2,o1) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o3,s_o4) == IEOR(s_o2,s_a) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o4,s_a) .and. &
IEOR(s_o2,s_o1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_(s_o2, s_o4, s_o3)%array(i_o2, i_o4, i_o3) =  &
    X_(s_o2, s_o4, s_o3)%array(i_o2, i_o4, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_o4, s_o3, s_o1)%array(i_o4, i_o3, i_o1) & 
  * h_(s_o1, s_o2)%array(i_o1, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x2(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x2(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x2(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o3, i_o3, s_o4, i_o4
integer :: s_o2, i_o2
! S2(i,k,m,a)  <-- 
! (   -1.00000000) D3(i,m,k,o3,o4,o2) X(o3,o4,o2,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o3,s_o4),s_o2) .and. &
IEOR(s_o3,s_o4) == IEOR(s_o2,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_o2, s_o4, s_o3, s_k, s_m, s_i)%array(i_o2, i_o4, i_o3, i_k, i_m, i_i) & 
  * X_(s_o2, s_o4, s_o3)%array(i_o2, i_o4, i_o3)
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

end subroutine g_sigma_ooov_ooov_no1_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x3(sa, ia, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x3(sa, ia, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x3(s_a, i_a, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_m, i_m, s_o3, i_o3
! X(o2,m,o3,a)  <-- 
! (    1.00000000)  T2(o1,o2,m,a) h(o3,o1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o2,s_m) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_m,s_a) .and. &
IEOR(s_o3,s_o1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_m, s_o2)%array(i_o3, i_m, i_o2) =  &
    X_(s_o3, s_m, s_o2)%array(i_o3, i_m, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_m, s_o2, s_o1)%array(i_m, i_o2, i_o1) & 
  * h_(s_o1, s_o3)%array(i_o1, i_o3)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x3(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x3(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x3(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_m, i_m
! S2(i,k,m,a)  <-- 
! (   -1.00000000) D2(i,o3,k,o2) X(o2,m,o3,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o2,s_m) == IEOR(s_o3,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * X_(s_o3, s_m, s_o2)%array(i_o3, i_m, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x4(sa, ia, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x4(sa, ia, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x4(s_a, i_a, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_o4, i_o4, s_o3, i_o3
! X(o2,o4,o3,a)  <-- 
! (    1.00000000)  T2(o2,o1,o4,a) h(o3,o1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o2,s_o4) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_o4,s_a) .and. &
IEOR(s_o3,s_o1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o4, s_o2)%array(i_o3, i_o4, i_o2) =  &
    X_(s_o3, s_o4, s_o2)%array(i_o3, i_o4, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_o4, s_o1, s_o2)%array(i_o4, i_o1, i_o2) & 
  * h_(s_o1, s_o3)%array(i_o1, i_o3)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x4(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x4(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x4(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o3, i_o3, s_o4, i_o4
integer :: s_o2, i_o2
! S2(i,k,m,a)  <-- 
! (   -1.00000000) D3(i,m,k,o3,o4,o2) X(o2,o4,o3,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o3,s_o4),s_o2) .and. &
IEOR(s_o2,s_o4) == IEOR(s_o3,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_o2, s_o4, s_o3, s_k, s_m, s_i)%array(i_o2, i_o4, i_o3, i_k, i_m, i_i) & 
  * X_(s_o3, s_o4, s_o2)%array(i_o3, i_o4, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x5(sa, ia, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x5(sa, ia, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x5(s_a, i_a, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o1, i_o1, s_m, i_m, s_o2, i_o2
! X(o3,m,o2,a)  <-- 
! (    1.00000000)  T2(o3,o1,m,a) h(o2,o1) 
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o3,s_m) == IEOR(s_o2,s_a) .and. & 
IEOR(s_o3,s_o1) == IEOR(s_m,s_a) .and. &
IEOR(s_o2,s_o1) == 0) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_(s_o2, s_m, s_o3)%array(i_o2, i_m, i_o3) =  &
    X_(s_o2, s_m, s_o3)%array(i_o2, i_m, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_m, s_o1, s_o3)%array(i_m, i_o1, i_o3) & 
  * h_(s_o1, s_o2)%array(i_o1, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x5(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x5(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x5(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_m, i_m
! S2(i,k,m,a)  <-- 
! (   -1.00000000) D2(i,o3,k,o2) X(o3,m,o2,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o3,s_m) == IEOR(s_o2,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * X_(s_o2, s_m, s_o3)%array(i_o2, i_m, i_o3)
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

end subroutine g_sigma_ooov_ooov_no1_x5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x6(sa, ia, sv1, iv1, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x6(sa, ia, sv1, iv1, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x6(s_a, i_a, s_v1, i_v1, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_o1, i_o1
! X(o3,o2,o1,a)  <-- 
! (    1.00000000)  T2(o3,o2,o1,v1) h(a,v1) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_o3,s_o2) == IEOR(s_o1,s_a) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_o1,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_o2, s_o3)%array(i_o1, i_o2, i_o3) =  &
    X_(s_o1, s_o2, s_o3)%array(i_o1, i_o2, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o2, s_o3)%array(i_o1, i_o2, i_o3) & 
  * h_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x6(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x6(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x6(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o2, i_o2, s_o1, i_o1
integer :: s_o3, i_o3
! S2(i,k,m,a)  <-- 
! (    1.00000000) D3(i,m,k,o2,o1,o3) X(o3,o2,o1,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o2,s_o1),s_o3) .and. &
IEOR(s_o3,s_o2) == IEOR(s_o1,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o3, s_o1, s_o2, s_k, s_m, s_i)%array(i_o3, i_o1, i_o2, i_k, i_m, i_i) & 
  * X_(s_o1, s_o2, s_o3)%array(i_o1, i_o2, i_o3)
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

end subroutine g_sigma_ooov_ooov_no1_x6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x7(sa, ia, sv1, iv1, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x7(sa, ia, sv1, iv1, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x7(s_a, i_a, s_v1, i_v1, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_m, i_m
! X(o1,o2,m,a)  <-- 
! (    1.00000000)  T2(o1,o2,m,v1) h(a,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_m,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_m,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o2, s_o1)%array(i_m, i_o2, i_o1) =  &
    X_(s_m, s_o2, s_o1)%array(i_m, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_m, s_o2, s_o1)%array(i_m, i_o2, i_o1) & 
  * h_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x7(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x7(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x7(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2, s_m, i_m
! S2(i,k,m,a)  <-- 
! (    1.00000000) D2(i,o1,k,o2) X(o1,o2,m,a) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_o2) .and. &
IEOR(s_o1,s_o2) == IEOR(s_m,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i) & 
  * X_(s_m, s_o2, s_o1)%array(i_m, i_o2, i_o1)
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

end subroutine g_sigma_ooov_ooov_no1_x7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y8(sc1, ic1, V2, Y0, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y0, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ooov_ooov_y8(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ooov_ooov_y8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y8(s_c1, i_c1, V2_, Y0_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y0_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o4, i_o4
! Y0(o1,o4)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o4) 
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o1,s_o4) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o4)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Y0_(s_o4, s_o1)%array(i_o4, i_o1) =  &
    Y0_(s_o4, s_o1)%array(i_o4, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o4, s_o1, s_c1)%array(i_o4, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x8(sa, ia, T2, Y0, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), Y0(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y0, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x8(sa, ia, av2_i, Yaa, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yaa)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x8(s_a, i_a, T2_, Y0_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y0_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o3, i_o3, s_o1, i_o1, s_o4, i_o4
! X(o2,o3,o4,a)  <-- 
! (    1.00000000)  T2(o2,o3,o1,a) Y0(o1,o4) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o2,s_o3) == IEOR(s_o4,s_a) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_o1,s_a) .and. &
IEOR(s_o1,s_o4) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o3, s_o2)%array(i_o4, i_o3, i_o2) =  &
    X_(s_o4, s_o3, s_o2)%array(i_o4, i_o3, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2) & 
  * Y0_(s_o4, s_o1)%array(i_o4, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x8(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x8(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x8(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o3, i_o3, s_o4, i_o4
integer :: s_o2, i_o2
! S2(i,k,m,a)  <-- 
! (    2.00000000) D3(i,m,k,o3,o4,o2) X(o2,o3,o4,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o3,s_o4),s_o2) .and. &
IEOR(s_o2,s_o3) == IEOR(s_o4,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 2.00000000d+00 & 
  * D3_(s_o2, s_o4, s_o3, s_k, s_m, s_i)%array(i_o2, i_o4, i_o3, i_k, i_m, i_i) & 
  * X_(s_o4, s_o3, s_o2)%array(i_o4, i_o3, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y9(sc1, ic1, V2, Y1, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y1(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y1, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ooov_ooov_y9(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ooov_ooov_y9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y9(s_c1, i_c1, V2_, Y1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o4, i_o4
! Y1(o1,o4)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o4) 
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o1,s_o4) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o4)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Y1_(s_o4, s_o1)%array(i_o4, i_o1) =  &
    Y1_(s_o4, s_o1)%array(i_o4, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o4, s_c1, s_o1)%array(i_o4, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x9(sa, ia, T2, Y1, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), Y1(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y1, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x9(sa, ia, av2_i, Yaa, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yaa)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x9(s_a, i_a, T2_, Y1_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y1_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o3, i_o3, s_o1, i_o1, s_o4, i_o4
! X(o2,o3,o4,a)  <-- 
! (    1.00000000)  T2(o2,o3,o1,a) Y1(o1,o4) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o2,s_o3) == IEOR(s_o4,s_a) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_o1,s_a) .and. &
IEOR(s_o1,s_o4) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o3, s_o2)%array(i_o4, i_o3, i_o2) =  &
    X_(s_o4, s_o3, s_o2)%array(i_o4, i_o3, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2) & 
  * Y1_(s_o4, s_o1)%array(i_o4, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x9(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x9(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x9(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o3, i_o3, s_o4, i_o4
integer :: s_o2, i_o2
! S2(i,k,m,a)  <-- 
! (   -1.00000000) D3(i,m,k,o3,o4,o2) X(o2,o3,o4,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o3,s_o4),s_o2) .and. &
IEOR(s_o2,s_o3) == IEOR(s_o4,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_o2, s_o4, s_o3, s_k, s_m, s_i)%array(i_o2, i_o4, i_o3, i_k, i_m, i_i) & 
  * X_(s_o4, s_o3, s_o2)%array(i_o4, i_o3, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y10(sc1, ic1, V2, Y2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y2, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ooov_ooov_y10(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ooov_ooov_y10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y10(s_c1, i_c1, V2_, Y2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y2_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_m, i_m, s_o1, i_o1
! Y2(m,o1)  <-- 
! (    1.00000000)  V2(c1,c1,m,o1) 
do s_m = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_m,s_o1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_m,s_o1)) then
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Y2_(s_o1, s_m)%array(i_o1, i_m) =  &
    Y2_(s_o1, s_m)%array(i_o1, i_m) &
  + 1.00000000d+00 & 
  * V2_(s_o1, s_m, s_c1)%array(i_o1, i_m, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x10(sa, ia, T2, Y2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), Y2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y2, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x10(sa, ia, av2_i, Yaa, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yaa)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x10(s_a, i_a, T2_, Y2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y2_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_o1, i_o1, s_m, i_m
! X(o3,o2,m,a)  <-- 
! (    1.00000000)  T2(o3,o2,o1,a) Y2(m,o1) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o3,s_o2) == IEOR(s_m,s_a) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_o1,s_a) .and. &
IEOR(s_m,s_o1) == 0) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o2, s_o3)%array(i_m, i_o2, i_o3) =  &
    X_(s_m, s_o2, s_o3)%array(i_m, i_o2, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o2, s_o3)%array(i_o1, i_o2, i_o3) & 
  * Y2_(s_o1, s_m)%array(i_o1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x10(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x10(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x10(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_m, i_m
! S2(i,k,m,a)  <-- 
! (    2.00000000) D2(i,o3,k,o2) X(o3,o2,m,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o3,s_o2) == IEOR(s_m,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 2.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * X_(s_m, s_o2, s_o3)%array(i_m, i_o2, i_o3)
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

end subroutine g_sigma_ooov_ooov_no1_x10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y11(sc1, ic1, V2, Y3, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y3(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y3, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ooov_ooov_y11(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ooov_ooov_y11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y11(s_c1, i_c1, V2_, Y3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y3_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_m, i_m, s_o1, i_o1
! Y3(m,o1)  <-- 
! (    1.00000000)  V2(c1,m,c1,o1) 
do s_m = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_m,s_o1) == 0 .and. & 
IEOR(s_c1,s_m) == IEOR(s_c1,s_o1)) then
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Y3_(s_o1, s_m)%array(i_o1, i_m) =  &
    Y3_(s_o1, s_m)%array(i_o1, i_m) &
  + 1.00000000d+00 & 
  * V2_(s_o1, s_c1, s_m)%array(i_o1, i_c1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x11(sa, ia, T2, Y3, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), Y3(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y3, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x11(sa, ia, av2_i, Yaa, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yaa)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x11(s_a, i_a, T2_, Y3_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y3_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_o1, i_o1, s_m, i_m
! X(o3,o2,m,a)  <-- 
! (    1.00000000)  T2(o3,o2,o1,a) Y3(m,o1) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o3,s_o2) == IEOR(s_m,s_a) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_o1,s_a) .and. &
IEOR(s_m,s_o1) == 0) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o2, s_o3)%array(i_m, i_o2, i_o3) =  &
    X_(s_m, s_o2, s_o3)%array(i_m, i_o2, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o2, s_o3)%array(i_o1, i_o2, i_o3) & 
  * Y3_(s_o1, s_m)%array(i_o1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x11(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x11(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x11(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_m, i_m
! S2(i,k,m,a)  <-- 
! (   -1.00000000) D2(i,o3,k,o2) X(o3,o2,m,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o3,s_o2) == IEOR(s_m,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * X_(s_m, s_o2, s_o3)%array(i_m, i_o2, i_o3)
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

end subroutine g_sigma_ooov_ooov_no1_x11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y12(sc1, ic1, V2, Y4, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y4(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y4, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ooov_ooov_y12(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ooov_ooov_y12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y12(s_c1, i_c1, V2_, Y4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y4_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y4(o1,o2)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y4_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y4_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_o1, s_c1)%array(i_o2, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x12(sa, ia, T2, Y4, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), Y4(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y4, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x12(sa, ia, av2_i, Yaa, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yaa)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x12(s_a, i_a, T2_, Y4_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y4_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3, s_o4, i_o4, s_o2, i_o2
! X(o3,o4,o2,a)  <-- 
! (    1.00000000)  T2(o1,o3,o4,a) Y4(o1,o2) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o3,s_o4) == IEOR(s_o2,s_a) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o4,s_a) .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_(s_o2, s_o4, s_o3)%array(i_o2, i_o4, i_o3) =  &
    X_(s_o2, s_o4, s_o3)%array(i_o2, i_o4, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_o4, s_o3, s_o1)%array(i_o4, i_o3, i_o1) & 
  * Y4_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x12(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x12(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x12(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o3, i_o3, s_o4, i_o4
integer :: s_o2, i_o2
! S2(i,k,m,a)  <-- 
! (   -2.00000000) D3(i,m,k,o3,o4,o2) X(o3,o4,o2,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o3,s_o4),s_o2) .and. &
IEOR(s_o3,s_o4) == IEOR(s_o2,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 2.00000000d+00 & 
  * D3_(s_o2, s_o4, s_o3, s_k, s_m, s_i)%array(i_o2, i_o4, i_o3, i_k, i_m, i_i) & 
  * X_(s_o2, s_o4, s_o3)%array(i_o2, i_o4, i_o3)
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

end subroutine g_sigma_ooov_ooov_no1_x12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y13(sc1, ic1, V2, Y5, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y5(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y5, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ooov_ooov_y13(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ooov_ooov_y13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y13(s_c1, i_c1, V2_, Y5_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y5_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y5(o1,o2)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y5_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y5_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_c1, s_o1)%array(i_o2, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x13(sa, ia, T2, Y5, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), Y5(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y5, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x13(sa, ia, av2_i, Yaa, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yaa)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x13(s_a, i_a, T2_, Y5_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y5_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3, s_o4, i_o4, s_o2, i_o2
! X(o3,o4,o2,a)  <-- 
! (    1.00000000)  T2(o1,o3,o4,a) Y5(o1,o2) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o3,s_o4) == IEOR(s_o2,s_a) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o4,s_a) .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_(s_o2, s_o4, s_o3)%array(i_o2, i_o4, i_o3) =  &
    X_(s_o2, s_o4, s_o3)%array(i_o2, i_o4, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_o4, s_o3, s_o1)%array(i_o4, i_o3, i_o1) & 
  * Y5_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x13(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x13(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x13(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o3, i_o3, s_o4, i_o4
integer :: s_o2, i_o2
! S2(i,k,m,a)  <-- 
! (    1.00000000) D3(i,m,k,o3,o4,o2) X(o3,o4,o2,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o3,s_o4),s_o2) .and. &
IEOR(s_o3,s_o4) == IEOR(s_o2,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o4, s_o3, s_k, s_m, s_i)%array(i_o2, i_o4, i_o3, i_k, i_m, i_i) & 
  * X_(s_o2, s_o4, s_o3)%array(i_o2, i_o4, i_o3)
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

end subroutine g_sigma_ooov_ooov_no1_x13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y14(sc1, ic1, V2, Y6, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y6(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y6, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ooov_ooov_y14(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ooov_ooov_y14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y14(s_c1, i_c1, V2_, Y6_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y6_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3
! Y6(o1,o3)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o3) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o1,s_o3) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o3)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
Y6_(s_o3, s_o1)%array(i_o3, i_o1) =  &
    Y6_(s_o3, s_o1)%array(i_o3, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o3, s_o1, s_c1)%array(i_o3, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x14(sa, ia, T2, Y6, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), Y6(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y6, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x14(sa, ia, av2_i, Yaa, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yaa)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x14(s_a, i_a, T2_, Y6_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y6_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_m, i_m, s_o3, i_o3
! X(o2,m,o3,a)  <-- 
! (    1.00000000)  T2(o1,o2,m,a) Y6(o1,o3) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o2,s_m) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_m,s_a) .and. &
IEOR(s_o1,s_o3) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_m, s_o2)%array(i_o3, i_m, i_o2) =  &
    X_(s_o3, s_m, s_o2)%array(i_o3, i_m, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_m, s_o2, s_o1)%array(i_m, i_o2, i_o1) & 
  * Y6_(s_o3, s_o1)%array(i_o3, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x14(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x14(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x14(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_m, i_m
! S2(i,k,m,a)  <-- 
! (   -2.00000000) D2(i,o3,k,o2) X(o2,m,o3,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o2,s_m) == IEOR(s_o3,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 2.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * X_(s_o3, s_m, s_o2)%array(i_o3, i_m, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y15(sc1, ic1, V2, Y7, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y7(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y7, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ooov_ooov_y15(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ooov_ooov_y15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y15(s_c1, i_c1, V2_, Y7_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y7_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3
! Y7(o1,o3)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o3) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o1,s_o3) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o3)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
Y7_(s_o3, s_o1)%array(i_o3, i_o1) =  &
    Y7_(s_o3, s_o1)%array(i_o3, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o3, s_c1, s_o1)%array(i_o3, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x15(sa, ia, T2, Y7, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), Y7(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y7, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x15(sa, ia, av2_i, Yaa, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yaa)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x15(s_a, i_a, T2_, Y7_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y7_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_m, i_m, s_o3, i_o3
! X(o2,m,o3,a)  <-- 
! (    1.00000000)  T2(o1,o2,m,a) Y7(o1,o3) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o2,s_m) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_m,s_a) .and. &
IEOR(s_o1,s_o3) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_m, s_o2)%array(i_o3, i_m, i_o2) =  &
    X_(s_o3, s_m, s_o2)%array(i_o3, i_m, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_m, s_o2, s_o1)%array(i_m, i_o2, i_o1) & 
  * Y7_(s_o3, s_o1)%array(i_o3, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x15(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x15(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x15(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_m, i_m
! S2(i,k,m,a)  <-- 
! (    1.00000000) D2(i,o3,k,o2) X(o2,m,o3,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o2,s_m) == IEOR(s_o3,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * X_(s_o3, s_m, s_o2)%array(i_o3, i_m, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y16(sc1, ic1, V2, Y8, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y8(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y8, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ooov_ooov_y16(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ooov_ooov_y16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y16(s_c1, i_c1, V2_, Y8_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y8_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3
! Y8(o1,o3)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o3) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o1,s_o3) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o3)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
Y8_(s_o3, s_o1)%array(i_o3, i_o1) =  &
    Y8_(s_o3, s_o1)%array(i_o3, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o3, s_o1, s_c1)%array(i_o3, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x16(sa, ia, T2, Y8, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), Y8(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y8, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x16(sa, ia, av2_i, Yaa, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yaa)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x16(s_a, i_a, T2_, Y8_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y8_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_o4, i_o4, s_o3, i_o3
! X(o2,o4,o3,a)  <-- 
! (    1.00000000)  T2(o2,o1,o4,a) Y8(o1,o3) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o2,s_o4) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_o4,s_a) .and. &
IEOR(s_o1,s_o3) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o4, s_o2)%array(i_o3, i_o4, i_o2) =  &
    X_(s_o3, s_o4, s_o2)%array(i_o3, i_o4, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_o4, s_o1, s_o2)%array(i_o4, i_o1, i_o2) & 
  * Y8_(s_o3, s_o1)%array(i_o3, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x16(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x16(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x16(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o3, i_o3, s_o4, i_o4
integer :: s_o2, i_o2
! S2(i,k,m,a)  <-- 
! (   -2.00000000) D3(i,m,k,o3,o4,o2) X(o2,o4,o3,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o3,s_o4),s_o2) .and. &
IEOR(s_o2,s_o4) == IEOR(s_o3,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 2.00000000d+00 & 
  * D3_(s_o2, s_o4, s_o3, s_k, s_m, s_i)%array(i_o2, i_o4, i_o3, i_k, i_m, i_i) & 
  * X_(s_o3, s_o4, s_o2)%array(i_o3, i_o4, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y17(sc1, ic1, V2, Y9, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y9(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y9, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ooov_ooov_y17(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ooov_ooov_y17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y17(s_c1, i_c1, V2_, Y9_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y9_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3
! Y9(o1,o3)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o3) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o1,s_o3) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o3)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
Y9_(s_o3, s_o1)%array(i_o3, i_o1) =  &
    Y9_(s_o3, s_o1)%array(i_o3, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o3, s_c1, s_o1)%array(i_o3, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x17(sa, ia, T2, Y9, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), Y9(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y9, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x17(sa, ia, av2_i, Yaa, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yaa)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x17(s_a, i_a, T2_, Y9_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y9_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_o4, i_o4, s_o3, i_o3
! X(o2,o4,o3,a)  <-- 
! (    1.00000000)  T2(o2,o1,o4,a) Y9(o1,o3) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o2,s_o4) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_o4,s_a) .and. &
IEOR(s_o1,s_o3) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o4, s_o2)%array(i_o3, i_o4, i_o2) =  &
    X_(s_o3, s_o4, s_o2)%array(i_o3, i_o4, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_o4, s_o1, s_o2)%array(i_o4, i_o1, i_o2) & 
  * Y9_(s_o3, s_o1)%array(i_o3, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x17(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x17(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x17(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o3, i_o3, s_o4, i_o4
integer :: s_o2, i_o2
! S2(i,k,m,a)  <-- 
! (    1.00000000) D3(i,m,k,o3,o4,o2) X(o2,o4,o3,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o3,s_o4),s_o2) .and. &
IEOR(s_o2,s_o4) == IEOR(s_o3,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o4, s_o3, s_k, s_m, s_i)%array(i_o2, i_o4, i_o3, i_k, i_m, i_i) & 
  * X_(s_o3, s_o4, s_o2)%array(i_o3, i_o4, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y18(sc1, ic1, V2, Y10, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y10(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y10, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ooov_ooov_y18(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ooov_ooov_y18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y18(s_c1, i_c1, V2_, Y10_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y10_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y10(o1,o2)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y10_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y10_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_o1, s_c1)%array(i_o2, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x18(sa, ia, T2, Y10, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), Y10(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y10, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x18(sa, ia, av2_i, Yaa, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yaa)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x18(s_a, i_a, T2_, Y10_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y10_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o1, i_o1, s_m, i_m, s_o2, i_o2
! X(o3,m,o2,a)  <-- 
! (    1.00000000)  T2(o3,o1,m,a) Y10(o1,o2) 
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o3,s_m) == IEOR(s_o2,s_a) .and. & 
IEOR(s_o3,s_o1) == IEOR(s_m,s_a) .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_(s_o2, s_m, s_o3)%array(i_o2, i_m, i_o3) =  &
    X_(s_o2, s_m, s_o3)%array(i_o2, i_m, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_m, s_o1, s_o3)%array(i_m, i_o1, i_o3) & 
  * Y10_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x18(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x18(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x18(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_m, i_m
! S2(i,k,m,a)  <-- 
! (   -2.00000000) D2(i,o3,k,o2) X(o3,m,o2,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o3,s_m) == IEOR(s_o2,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 2.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * X_(s_o2, s_m, s_o3)%array(i_o2, i_m, i_o3)
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

end subroutine g_sigma_ooov_ooov_no1_x18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y19(sc1, ic1, V2, Y11, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y11(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y11, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ooov_ooov_y19(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ooov_ooov_y19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y19(s_c1, i_c1, V2_, Y11_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y11_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y11(o1,o2)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y11_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y11_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_c1, s_o1)%array(i_o2, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x19(sa, ia, T2, Y11, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), Y11(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y11, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x19(sa, ia, av2_i, Yaa, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yaa)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x19(s_a, i_a, T2_, Y11_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y11_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o1, i_o1, s_m, i_m, s_o2, i_o2
! X(o3,m,o2,a)  <-- 
! (    1.00000000)  T2(o3,o1,m,a) Y11(o1,o2) 
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o3,s_m) == IEOR(s_o2,s_a) .and. & 
IEOR(s_o3,s_o1) == IEOR(s_m,s_a) .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_(s_o2, s_m, s_o3)%array(i_o2, i_m, i_o3) =  &
    X_(s_o2, s_m, s_o3)%array(i_o2, i_m, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_m, s_o1, s_o3)%array(i_m, i_o1, i_o3) & 
  * Y11_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x19(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x19(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x19(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_m, i_m
! S2(i,k,m,a)  <-- 
! (    1.00000000) D2(i,o3,k,o2) X(o3,m,o2,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o3,s_m) == IEOR(s_o2,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * X_(s_o2, s_m, s_o3)%array(i_o2, i_m, i_o3)
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

end subroutine g_sigma_ooov_ooov_no1_x19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y20(sc1, ic1, V2, Y12, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y12(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y12, nir, nsym, psym) ! -> Yvv (allocate) 
call g_sigma_ooov_ooov_y20(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ooov_ooov_y20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y20(s_c1, i_c1, V2_, Y12_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y12_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_v1, i_v1
! Y12(a,v1)  <-- 
! (    1.00000000)  V2(c1,c1,a,v1) 
do s_a = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_v1)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y12_(s_v1, s_a)%array(i_v1, i_a) =  &
    Y12_(s_v1, s_a)%array(i_v1, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_a, s_c1)%array(i_v1, i_a, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x20(sa, ia, sv1, iv1, T2, Y12, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y12(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y12, nir, nsym, psym) ! -> Yvv (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x20(sa, ia, sv1, iv1, av2_i, Yvv, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x20(s_a, i_a, s_v1, i_v1, T2_, Y12_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y12_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_o3, i_o3
! X(o1,o2,o3,a)  <-- 
! (    1.00000000)  T2(o1,o2,o3,v1) Y12(a,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_o3,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1) =  &
    X_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1) & 
  * Y12_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x20(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x20(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x20(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o2, i_o2, s_o3, i_o3
integer :: s_o1, i_o1
! S2(i,k,m,a)  <-- 
! (    2.00000000) D3(i,m,k,o2,o3,o1) X(o1,o2,o3,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o2,s_o3),s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_o3,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 2.00000000d+00 & 
  * D3_(s_o1, s_o3, s_o2, s_k, s_m, s_i)%array(i_o1, i_o3, i_o2, i_k, i_m, i_i) & 
  * X_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1)
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

end subroutine g_sigma_ooov_ooov_no1_x20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y21(sc1, ic1, V2, Y13, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y13(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y13, nir, nsym, psym) ! -> Yvv (allocate) 
call g_sigma_ooov_ooov_y21(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ooov_ooov_y21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y21(s_c1, i_c1, V2_, Y13_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y13_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_v1, i_v1
! Y13(a,v1)  <-- 
! (    1.00000000)  V2(c1,a,c1,v1) 
do s_a = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_v1)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y13_(s_v1, s_a)%array(i_v1, i_a) =  &
    Y13_(s_v1, s_a)%array(i_v1, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_a)%array(i_v1, i_c1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x21(sa, ia, sv1, iv1, T2, Y13, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y13(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y13, nir, nsym, psym) ! -> Yvv (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x21(sa, ia, sv1, iv1, av2_i, Yvv, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x21(s_a, i_a, s_v1, i_v1, T2_, Y13_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y13_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_o3, i_o3
! X(o1,o2,o3,a)  <-- 
! (    1.00000000)  T2(o1,o2,o3,v1) Y13(a,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_o3,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1) =  &
    X_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1) & 
  * Y13_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x21(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x21(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x21(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o2, i_o2, s_o3, i_o3
integer :: s_o1, i_o1
! S2(i,k,m,a)  <-- 
! (   -1.00000000) D3(i,m,k,o2,o3,o1) X(o1,o2,o3,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o2,s_o3),s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_o3,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_o1, s_o3, s_o2, s_k, s_m, s_i)%array(i_o1, i_o3, i_o2, i_k, i_m, i_i) & 
  * X_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1)
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

end subroutine g_sigma_ooov_ooov_no1_x21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y22(sc1, ic1, V2, Y14, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y14(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y14, nir, nsym, psym) ! -> Yvv (allocate) 
call g_sigma_ooov_ooov_y22(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ooov_ooov_y22



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y22(s_c1, i_c1, V2_, Y14_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y14_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_v1, i_v1
! Y14(a,v1)  <-- 
! (    1.00000000)  V2(c1,c1,a,v1) 
do s_a = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_v1)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y14_(s_v1, s_a)%array(i_v1, i_a) =  &
    Y14_(s_v1, s_a)%array(i_v1, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_a, s_c1)%array(i_v1, i_a, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y22



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x22(sa, ia, sv1, iv1, T2, Y14, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y14(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y14, nir, nsym, psym) ! -> Yvv (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x22(sa, ia, sv1, iv1, av2_i, Yvv, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x22



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x22(s_a, i_a, s_v1, i_v1, T2_, Y14_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y14_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_m, i_m
! X(o2,o1,m,a)  <-- 
! (    1.00000000)  T2(o2,o1,m,v1) Y14(a,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_m,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_m,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) =  &
    X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) & 
  * Y14_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x22



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x22(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x22(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x22



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x22(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o1, i_o1, s_m, i_m
! S2(i,k,m,a)  <-- 
! (    2.00000000) D2(i,o2,k,o1) X(o2,o1,m,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o2,s_o1) == IEOR(s_m,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 2.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
  * X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x22



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_y23(sc1, ic1, V2, Y15, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y15(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y15, nir, nsym, psym) ! -> Yvv (allocate) 
call g_sigma_ooov_ooov_y23(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ooov_ooov_y23



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_y23(s_c1, i_c1, V2_, Y15_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y15_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_v1, i_v1
! Y15(a,v1)  <-- 
! (    1.00000000)  V2(c1,a,c1,v1) 
do s_a = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_v1)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y15_(s_v1, s_a)%array(i_v1, i_a) =  &
    Y15_(s_v1, s_a)%array(i_v1, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_a)%array(i_v1, i_c1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_y23



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x23(sa, ia, sv1, iv1, T2, Y15, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y15(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y15, nir, nsym, psym) ! -> Yvv (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x23(sa, ia, sv1, iv1, av2_i, Yvv, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x23



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x23(s_a, i_a, s_v1, i_v1, T2_, Y15_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y15_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_m, i_m
! X(o2,o1,m,a)  <-- 
! (    1.00000000)  T2(o2,o1,m,v1) Y15(a,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_m,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_m,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) =  &
    X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) & 
  * Y15_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x23



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x23(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x23(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x23



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x23(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o1, i_o1, s_m, i_m
! S2(i,k,m,a)  <-- 
! (   -1.00000000) D2(i,o2,k,o1) X(o2,o1,m,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o2,s_o1) == IEOR(s_m,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
  * X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x23



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x24(so3, io3, so6, io6, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: so3, io3, so6, io6
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so6, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = so3

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x24(so3, io3, so6, io6, h2_i, xaaaaa, d4_ij, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x24



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x24(s_o3, i_o3, s_o6, i_o6, V2_, X_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o6, s_o6
integer, intent(in) :: i_o3, s_o3
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o4, i_o4, s_o5, i_o5
integer :: s_o2, i_o2, s_o1, i_o1
! X(o3,i,m,k,o4,o1)  <-- 
! (    1.00000000)  D4(o6,o3,i,m,k,o4,o5,o2) V2(o6,o1,o2,o5) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o4 = 0, nir-1
do s_o5 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(IEOR(s_o3,s_i),s_m) == IEOR(IEOR(s_k,s_o4),s_o1) .and. & 
IEOR(IEOR(s_o6,s_o3),IEOR(s_i,s_m)) == IEOR(IEOR(s_k,s_o4),IEOR(s_o5,s_o2)) .and. &
IEOR(s_o6,s_o1) == IEOR(s_o2,s_o5)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_o4, s_k, s_m, s_i)%array(i_o1, i_o4, i_k, i_m, i_i) =  &
    X_(s_o1, s_o4, s_k, s_m, s_i)%array(i_o1, i_o4, i_k, i_m, i_i) &
  + 1.00000000d+00 & 
  * D4_(s_o2, s_o5, s_o4, s_k, s_m, s_i)%array(i_o2, i_o5, i_o4, i_k, i_m, i_i) & 
  * V2_(s_o5, s_o2, s_o1)%array(i_o5, i_o2, i_o1)
end do ! Irrep Loop
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
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x24



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x24(sa, ia, so3, io3, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, so3, io3
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = so3

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x24(sa, ia, so3, io3, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x24



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x24(s_a, i_a, s_o3, i_o3, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o3, s_o3
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o4, i_o4, s_o1, i_o1, s_i, i_i, s_m, i_m, s_k, i_k
! S2(i,k,m,a)  <-- 
! (    1.00000000) T2(o3,o4,o1,a) X(o3,i,m,k,o4,o1) 
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_o3,s_o4) == IEOR(s_o1,s_a) .and. &
IEOR(IEOR(s_o3,s_i),s_m) == IEOR(IEOR(s_k,s_o4),s_o1)) then
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o4, s_o3)%array(i_o1, i_o4, i_o3) & 
  * X_(s_o1, s_o4, s_k, s_m, s_i)%array(i_o1, i_o4, i_k, i_m, i_i)
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

end subroutine g_sigma_ooov_ooov_no1_x24



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x25(sm, im, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sm, im
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sm, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sm

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x25(sm, im, h2_i, xaaaaa, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x25



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x25(s_m, i_m, V2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_m, s_m
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o4, i_o4, s_k, i_k, s_o3, i_o3, s_o5, i_o5
integer :: s_o2, i_o2, s_o1, i_o1
! X(i,o4,k,o3,m,o1)  <-- 
! (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(m,o1,o2,o5) 
do s_i = 0, nir-1
do s_o4 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o5 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_m),s_o1) .and. & 
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_o5),s_o2) .and. &
IEOR(s_m,s_o1) == IEOR(s_o2,s_o5)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_o3, s_k, s_o4, s_i)%array(i_o1, i_o3, i_k, i_o4, i_i) =  &
    X_(s_o1, s_o3, s_k, s_o4, s_i)%array(i_o1, i_o3, i_k, i_o4, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o5, s_o3, s_k, s_o4, s_i)%array(i_o2, i_o5, i_o3, i_k, i_o4, i_i) & 
  * V2_(s_o5, s_o2, s_o1)%array(i_o5, i_o2, i_o1)
end do ! Irrep Loop
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
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x25



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x25(sa, ia, sm, im, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sm, im
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sm

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x25(sa, ia, sm, im, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x25



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x25(s_a, i_a, s_m, i_m, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_m, s_m
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o4, i_o4, s_o3, i_o3, s_o1, i_o1, s_i, i_i, s_k, i_k
! S2(i,k,m,a)  <-- 
! (    1.00000000) T2(o4,o3,o1,a) X(i,o4,k,o3,m,o1) 
do s_o4 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_o4,s_o3) == IEOR(s_o1,s_a) .and. &
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_m),s_o1)) then
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o3, s_o4)%array(i_o1, i_o3, i_o4) & 
  * X_(s_o1, s_o3, s_k, s_o4, s_i)%array(i_o1, i_o3, i_k, i_o4, i_i)
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

end subroutine g_sigma_ooov_ooov_no1_x25



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x26(sm, im, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sm, im
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sm, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sm

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x26(sm, im, h2_i, xaaaaa, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x26



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x26(s_m, i_m, V2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_m, s_m
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o4, i_o4, s_k, i_k, s_o3, i_o3, s_o5, i_o5
integer :: s_o2, i_o2, s_o1, i_o1
! X(i,k,o3,o2,m,o1)  <-- 
! (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(m,o4,o1,o5) 
do s_i = 0, nir-1
do s_o4 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o5 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(IEOR(s_i,s_k),s_o3) == IEOR(IEOR(s_o2,s_m),s_o1) .and. & 
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_o5),s_o2) .and. &
IEOR(s_m,s_o4) == IEOR(s_o1,s_o5)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_o2, s_o3, s_k, s_i)%array(i_o1, i_o2, i_o3, i_k, i_i) =  &
    X_(s_o1, s_o2, s_o3, s_k, s_i)%array(i_o1, i_o2, i_o3, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o5, s_o3, s_k, s_o4, s_i)%array(i_o2, i_o5, i_o3, i_k, i_o4, i_i) & 
  * V2_(s_o5, s_o1, s_o4)%array(i_o5, i_o1, i_o4)
end do ! Irrep Loop
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
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x26



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x26(sa, ia, sm, im, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sm, im
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sm

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x26(sa, ia, sm, im, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x26



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x26(s_a, i_a, s_m, i_m, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_m, s_m
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o3, i_o3, s_o1, i_o1, s_i, i_i, s_k, i_k
! S2(i,k,m,a)  <-- 
! (    1.00000000) T2(o2,o3,o1,a) X(i,k,o3,o2,m,o1) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_o1,s_a) .and. &
IEOR(IEOR(s_i,s_k),s_o3) == IEOR(IEOR(s_o2,s_m),s_o1)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2) & 
  * X_(s_o1, s_o2, s_o3, s_k, s_i)%array(i_o1, i_o2, i_o3, i_k, i_i)
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

end subroutine g_sigma_ooov_ooov_no1_x26



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x27(so3, io3, so6, io6, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: so3, io3, so6, io6
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so3, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = so6

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x27(so3, io3, so6, io6, h2_i, xaaaaa, d4_ij, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x27



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x27(s_o3, i_o3, s_o6, i_o6, V2_, X_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o3, s_o3
integer, intent(in) :: i_o6, s_o6
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_m, i_m, s_i, i_i, s_o4, i_o4, s_k, i_k, s_o2, i_o2
integer :: s_o5, i_o5, s_o1, i_o1
! X(o6,m,i,o4,k,o1)  <-- 
! (    1.00000000)  D4(o3,o6,m,i,o4,k,o2,o5) V2(o3,o1,o2,o5) 
do s_m = 0, nir-1
do s_i = 0, nir-1
do s_o4 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o5 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(IEOR(s_o6,s_m),s_i) == IEOR(IEOR(s_o4,s_k),s_o1) .and. & 
IEOR(IEOR(s_o3,s_o6),IEOR(s_m,s_i)) == IEOR(IEOR(s_o4,s_k),IEOR(s_o2,s_o5)) .and. &
IEOR(s_o3,s_o1) == IEOR(s_o2,s_o5)) then
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_k, s_o4, s_i, s_m)%array(i_o1, i_k, i_o4, i_i, i_m) =  &
    X_(s_o1, s_k, s_o4, s_i, s_m)%array(i_o1, i_k, i_o4, i_i, i_m) &
  + 1.00000000d+00 & 
  * D4_(s_o5, s_o2, s_k, s_o4, s_i, s_m)%array(i_o5, i_o2, i_k, i_o4, i_i, i_m) & 
  * V2_(s_o5, s_o2, s_o1)%array(i_o5, i_o2, i_o1)
end do ! Irrep Loop
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
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x27



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x27(sa, ia, so6, io6, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, so6, io6
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = so6

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x27(sa, ia, so6, io6, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x27



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x27(s_a, i_a, s_o6, i_o6, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o6, s_o6
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o4, i_o4, s_m, i_m, s_i, i_i, s_k, i_k
! S2(i,k,m,a)  <-- 
! (   -1.00000000) T2(o1,o4,o6,a) X(o6,m,i,o4,k,o1) 
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_m = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_o1,s_o4) == IEOR(s_o6,s_a) .and. &
IEOR(IEOR(s_o6,s_m),s_i) == IEOR(IEOR(s_o4,s_k),s_o1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * T2_(s_o6, s_o4, s_o1)%array(i_o6, i_o4, i_o1) & 
  * X_(s_o1, s_k, s_o4, s_i, s_m)%array(i_o1, i_k, i_o4, i_i, i_m)
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

end subroutine g_sigma_ooov_ooov_no1_x27



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x28(so1, io1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: so1, io1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = so1

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x28(so1, io1, h2_i, xaaa, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x28



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x28(s_o1, i_o1, V2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o4, i_o4, s_k, i_k, s_o3, i_o3, s_o5, i_o5
integer :: s_o2, i_o2
! X(i,k,o3,o1)  <-- 
! (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(o1,o4,o2,o5) 
do s_i = 0, nir-1
do s_o4 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o5 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_o3,s_o1) .and. & 
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_o5),s_o2) .and. &
IEOR(s_o1,s_o4) == IEOR(s_o2,s_o5)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_(s_o3, s_k, s_i)%array(i_o3, i_k, i_i) =  &
    X_(s_o3, s_k, s_i)%array(i_o3, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o5, s_o3, s_k, s_o4, s_i)%array(i_o2, i_o5, i_o3, i_k, i_o4, i_i) & 
  * V2_(s_o5, s_o2, s_o4)%array(i_o5, i_o2, i_o4)
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

end subroutine g_sigma_ooov_ooov_no0_x28



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x28(sa, ia, so1, io1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, so1, io1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = so1

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x28(sa, ia, so1, io1, av2_i, xaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x28



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x28(s_a, i_a, s_o1, i_o1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_m, i_m, s_i, i_i, s_k, i_k
! S2(i,k,m,a)  <-- 
! (   -1.00000000) T2(o1,o3,m,a) X(i,k,o3,o1) 
do s_o3 = 0, nir-1
do s_m = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_m,s_a) .and. &
IEOR(s_i,s_k) == IEOR(s_o3,s_o1)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * T2_(s_m, s_o3, s_o1)%array(i_m, i_o3, i_o1) & 
  * X_(s_o3, s_k, s_i)%array(i_o3, i_k, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no1_x28



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x29(sm, im, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sm, im
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sm, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sm

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x29(sm, im, h2_i, xaaaaa, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x29



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x29(s_m, i_m, V2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_m, s_m
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o4, i_o4, s_k, i_k, s_o3, i_o3, s_o5, i_o5
integer :: s_o2, i_o2, s_o1, i_o1
! X(i,k,o3,o5,m,o1)  <-- 
! (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(m,o4,o1,o2) 
do s_i = 0, nir-1
do s_o4 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o5 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(IEOR(s_i,s_k),s_o3) == IEOR(IEOR(s_o5,s_m),s_o1) .and. & 
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_o5),s_o2) .and. &
IEOR(s_m,s_o4) == IEOR(s_o1,s_o2)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_o5, s_o3, s_k, s_i)%array(i_o1, i_o5, i_o3, i_k, i_i) =  &
    X_(s_o1, s_o5, s_o3, s_k, s_i)%array(i_o1, i_o5, i_o3, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o5, s_o3, s_k, s_o4, s_i)%array(i_o2, i_o5, i_o3, i_k, i_o4, i_i) & 
  * V2_(s_o2, s_o1, s_o4)%array(i_o2, i_o1, i_o4)
end do ! Irrep Loop
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
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x29



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x29(sa, ia, sm, im, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sm, im
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sm

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x29(sa, ia, sm, im, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x29



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x29(s_a, i_a, s_m, i_m, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_m, s_m
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3, s_o5, i_o5, s_i, i_i, s_k, i_k
! S2(i,k,m,a)  <-- 
! (   -1.00000000) T2(o1,o3,o5,a) X(i,k,o3,o5,m,o1) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o5 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o5,s_a) .and. &
IEOR(IEOR(s_i,s_k),s_o3) == IEOR(IEOR(s_o5,s_m),s_o1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * T2_(s_o5, s_o3, s_o1)%array(i_o5, i_o3, i_o1) & 
  * X_(s_o1, s_o5, s_o3, s_k, s_i)%array(i_o1, i_o5, i_o3, i_k, i_i)
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

end subroutine g_sigma_ooov_ooov_no1_x29



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x30(sk, ik, so4, io4, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sk, ik, so4, io4
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so4, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sk

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x30(sk, ik, so4, io4, h2_i, xaaaaa, d4_ij, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x30



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x30(s_k, i_k, s_o4, i_o4, V2_, X_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o4, s_o4
integer, intent(in) :: i_k, s_k
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_m, i_m, s_i, i_i, s_o2, i_o2, s_o5, i_o5, s_o3, i_o3
integer :: s_o6, i_o6, s_o1, i_o1
! X(k,m,i,o3,o6,o1)  <-- 
! (    1.00000000)  D4(o4,k,m,i,o2,o5,o3,o6) V2(o4,o1,o2,o5) 
do s_m = 0, nir-1
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_o5 = 0, nir-1
do s_o3 = 0, nir-1
do s_o6 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(IEOR(s_k,s_m),s_i) == IEOR(IEOR(s_o3,s_o6),s_o1) .and. & 
IEOR(IEOR(s_o4,s_k),IEOR(s_m,s_i)) == IEOR(IEOR(s_o2,s_o5),IEOR(s_o3,s_o6)) .and. &
IEOR(s_o4,s_o1) == IEOR(s_o2,s_o5)) then
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o6 = psym(I_BEGIN, I_O, s_o6), psym(I_END, I_O, s_o6)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_o6, s_o3, s_i, s_m)%array(i_o1, i_o6, i_o3, i_i, i_m) =  &
    X_(s_o1, s_o6, s_o3, s_i, s_m)%array(i_o1, i_o6, i_o3, i_i, i_m) &
  + 1.00000000d+00 & 
  * D4_(s_o6, s_o3, s_o5, s_o2, s_i, s_m)%array(i_o6, i_o3, i_o5, i_o2, i_i, i_m) & 
  * V2_(s_o5, s_o2, s_o1)%array(i_o5, i_o2, i_o1)
end do ! Irrep Loop
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
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x30



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x30(sa, ia, sk, ik, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sk, ik
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sk

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x30(sa, ia, sk, ik, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x30



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x30(s_a, i_a, s_k, i_k, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_k, s_k
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o1, i_o1, s_o6, i_o6, s_m, i_m, s_i, i_i
! S2(i,k,m,a)  <-- 
! (   -1.00000000) T2(o3,o1,o6,a) X(k,m,i,o3,o6,o1) 
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_o6 = 0, nir-1
do s_m = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_o3,s_o1) == IEOR(s_o6,s_a) .and. &
IEOR(IEOR(s_k,s_m),s_i) == IEOR(IEOR(s_o3,s_o6),s_o1)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o6 = psym(I_BEGIN, I_O, s_o6), psym(I_END, I_O, s_o6)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * T2_(s_o6, s_o1, s_o3)%array(i_o6, i_o1, i_o3) & 
  * X_(s_o1, s_o6, s_o3, s_i, s_m)%array(i_o1, i_o6, i_o3, i_i, i_m)
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

end subroutine g_sigma_ooov_ooov_no1_x30



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x31(so1, io1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: so1, io1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = so1

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x31(so1, io1, h2_i, xaaa, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x31



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x31(s_o1, i_o1, V2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o4, i_o4, s_k, i_k, s_o3, i_o3, s_o5, i_o5
integer :: s_o2, i_o2
! X(i,o4,k,o1)  <-- 
! (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(o1,o3,o2,o5) 
do s_i = 0, nir-1
do s_o4 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o5 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_o4) == IEOR(s_k,s_o1) .and. & 
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_o5),s_o2) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o5)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_(s_k, s_o4, s_i)%array(i_k, i_o4, i_i) =  &
    X_(s_k, s_o4, s_i)%array(i_k, i_o4, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o5, s_o3, s_k, s_o4, s_i)%array(i_o2, i_o5, i_o3, i_k, i_o4, i_i) & 
  * V2_(s_o5, s_o2, s_o3)%array(i_o5, i_o2, i_o3)
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

end subroutine g_sigma_ooov_ooov_no0_x31



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x31(sa, ia, so1, io1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, so1, io1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = so1

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x31(sa, ia, so1, io1, av2_i, xaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x31



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x31(s_a, i_a, s_o1, i_o1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o4, i_o4, s_m, i_m, s_i, i_i, s_k, i_k
! S2(i,k,m,a)  <-- 
! (   -1.00000000) T2(o4,o1,m,a) X(i,o4,k,o1) 
do s_o4 = 0, nir-1
do s_m = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_o4,s_o1) == IEOR(s_m,s_a) .and. &
IEOR(s_i,s_o4) == IEOR(s_k,s_o1)) then
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * T2_(s_m, s_o1, s_o4)%array(i_m, i_o1, i_o4) & 
  * X_(s_k, s_o4, s_i)%array(i_k, i_o4, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no1_x31



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x32(sm, im, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sm, im
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sm, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sm

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x32(sm, im, h2_i, xaaaaa, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x32



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x32(s_m, i_m, V2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_m, s_m
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o4, i_o4, s_k, i_k, s_o3, i_o3, s_o5, i_o5
integer :: s_o2, i_o2, s_o1, i_o1
! X(i,k,o5,o2,m,o1)  <-- 
! (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(m,o4,o1,o3) 
do s_i = 0, nir-1
do s_o4 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o5 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(IEOR(s_i,s_k),s_o5) == IEOR(IEOR(s_o2,s_m),s_o1) .and. & 
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_o5),s_o2) .and. &
IEOR(s_m,s_o4) == IEOR(s_o1,s_o3)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_o2, s_o5, s_k, s_i)%array(i_o1, i_o2, i_o5, i_k, i_i) =  &
    X_(s_o1, s_o2, s_o5, s_k, s_i)%array(i_o1, i_o2, i_o5, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o5, s_o3, s_k, s_o4, s_i)%array(i_o2, i_o5, i_o3, i_k, i_o4, i_i) & 
  * V2_(s_o3, s_o1, s_o4)%array(i_o3, i_o1, i_o4)
end do ! Irrep Loop
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
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x32



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x32(sa, ia, sm, im, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sm, im
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sm

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x32(sa, ia, sm, im, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x32



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x32(s_a, i_a, s_m, i_m, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_m, s_m
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_o5, i_o5, s_i, i_i, s_k, i_k
! S2(i,k,m,a)  <-- 
! (   -1.00000000) T2(o2,o1,o5,a) X(i,k,o5,o2,m,o1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o5 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_o5,s_a) .and. &
IEOR(IEOR(s_i,s_k),s_o5) == IEOR(IEOR(s_o2,s_m),s_o1)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * T2_(s_o5, s_o1, s_o2)%array(i_o5, i_o1, i_o2) & 
  * X_(s_o1, s_o2, s_o5, s_k, s_i)%array(i_o1, i_o2, i_o5, i_k, i_i)
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

end subroutine g_sigma_ooov_ooov_no1_x32



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x33(sa, ia, so1, io1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, so1, io1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x33(sa, ia, so1, io1, av2_i, h2_i, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x33



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x33(s_a, i_a, s_o1, i_o1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o5, i_o5, s_o4, i_o4, s_o3, i_o3
! X(o5,o4,o3,a)  <-- 
! (    1.00000000)  T2(o2,o1,o5,a) V2(o1,o4,o2,o3) 
do s_o2 = 0, nir-1
do s_o5 = 0, nir-1
do s_o4 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o5,s_o4) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_o5,s_a) .and. &
IEOR(s_o1,s_o4) == IEOR(s_o2,s_o3)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o4, s_o5)%array(i_o3, i_o4, i_o5) =  &
    X_(s_o3, s_o4, s_o5)%array(i_o3, i_o4, i_o5) &
  + 1.00000000d+00 & 
  * T2_(s_o5, s_o1, s_o2)%array(i_o5, i_o1, i_o2) & 
  * V2_(s_o3, s_o2, s_o4)%array(i_o3, i_o2, i_o4)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x33



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x33(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x33(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x33



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x33(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o4, i_o4, s_o5, i_o5
integer :: s_o3, i_o3
! S2(i,k,m,a)  <-- 
! (   -1.00000000) D3(i,m,k,o4,o5,o3) X(o5,o4,o3,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o4 = 0, nir-1
do s_o5 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o4,s_o5),s_o3) .and. &
IEOR(s_o5,s_o4) == IEOR(s_o3,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_o3, s_o5, s_o4, s_k, s_m, s_i)%array(i_o3, i_o5, i_o4, i_k, i_m, i_i) & 
  * X_(s_o3, s_o4, s_o5)%array(i_o3, i_o4, i_o5)
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

end subroutine g_sigma_ooov_ooov_no1_x33



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x34(sa, ia, so1, io1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, so1, io1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x34(sa, ia, so1, io1, av2_i, h2_i, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x34



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x34(s_a, i_a, s_o1, i_o1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_m, i_m, s_o3, i_o3, s_o4, i_o4
! X(m,o3,o4,a)  <-- 
! (    1.00000000)  T2(o2,o1,m,a) V2(o1,o3,o2,o4) 
do s_o2 = 0, nir-1
do s_m = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_m,s_o3) == IEOR(s_o4,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_m,s_a) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o3, s_m)%array(i_o4, i_o3, i_m) =  &
    X_(s_o4, s_o3, s_m)%array(i_o4, i_o3, i_m) &
  + 1.00000000d+00 & 
  * T2_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) & 
  * V2_(s_o4, s_o2, s_o3)%array(i_o4, i_o2, i_o3)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x34



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x34(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x34(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x34



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x34(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o4, i_o4, s_k, i_k, s_o3, i_o3, s_m, i_m
! S2(i,k,m,a)  <-- 
! (   -1.00000000) D2(i,o4,k,o3) X(m,o3,o4,a) 
do s_i = 0, nir-1
do s_o4 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o4) == IEOR(s_k,s_o3) .and. &
IEOR(s_m,s_o3) == IEOR(s_o4,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_o3, s_k, s_o4, s_i)%array(i_o3, i_k, i_o4, i_i) & 
  * X_(s_o4, s_o3, s_m)%array(i_o4, i_o3, i_m)
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

end subroutine g_sigma_ooov_ooov_no1_x34



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x35(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x35(sa, ia, sv1, iv1, av2_i, h2_i, xaaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x35



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x35(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o3, i_o3, s_o5, i_o5, s_o1, i_o1, s_o4, i_o4
! X(o2,o3,o5,o1,o4,a)  <-- 
! (    1.00000000)  T2(o2,o3,o5,v1) V2(a,v1,o1,o4) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o5 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(IEOR(s_o2,s_o3),s_o5) == IEOR(IEOR(s_o1,s_o4),s_a) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_o5,s_v1) .and. &
IEOR(s_a,s_v1) == IEOR(s_o1,s_o4)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o1, s_o5, s_o3, s_o2)%array(i_o4, i_o1, i_o5, i_o3, i_o2) =  &
    X_(s_o4, s_o1, s_o5, s_o3, s_o2)%array(i_o4, i_o1, i_o5, i_o3, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_o5, s_o3, s_o2)%array(i_o5, i_o3, i_o2) & 
  * V2_(s_o4, s_o1, s_v1)%array(i_o4, i_o1, i_v1)
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

end subroutine g_sigma_ooov_ooov_no0_x35



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x35(sa, ia, si, ii, sm, im, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, si, ii, sm, im
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x35(sa, ia, si, ii, sm, im, xaaaaa, av2_i2, d4_ij, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x35



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x35(s_a, i_a, s_i, i_i, s_m, i_m, X_, S2_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_m, s_m
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_o3, i_o3, s_o4, i_o4, s_o1, i_o1, s_o5, i_o5
integer :: s_o2, i_o2
! S2(i,k,m,a)  <-- 
! (    1.00000000) D4(i,m,k,o3,o4,o1,o5,o2) X(o2,o3,o5,o1,o4,a) 
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
do s_o5 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),IEOR(s_k,s_o3)) == IEOR(IEOR(s_o4,s_o1),IEOR(s_o5,s_o2)) .and. &
IEOR(IEOR(s_o2,s_o3),s_o5) == IEOR(IEOR(s_o1,s_o4),s_a)) then
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D4_(s_o2, s_o5, s_o1, s_o4, s_o3, s_k)%array(i_o2, i_o5, i_o1, i_o4, i_o3, i_k) & 
  * X_(s_o4, s_o1, s_o5, s_o3, s_o2)%array(i_o4, i_o1, i_o5, i_o3, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x35



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x36(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x36(sa, ia, sv1, iv1, av2_i, h2_i, xaaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x36



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x36(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_m, i_m, s_o1, i_o1, s_o4, i_o4
! X(o3,o2,m,o1,o4,a)  <-- 
! (    1.00000000)  T2(o3,o2,m,v1) V2(a,v1,o1,o4) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(IEOR(s_o3,s_o2),s_m) == IEOR(IEOR(s_o1,s_o4),s_a) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_m,s_v1) .and. &
IEOR(s_a,s_v1) == IEOR(s_o1,s_o4)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o1, s_m, s_o2, s_o3)%array(i_o4, i_o1, i_m, i_o2, i_o3) =  &
    X_(s_o4, s_o1, s_m, s_o2, s_o3)%array(i_o4, i_o1, i_m, i_o2, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_m, s_o2, s_o3)%array(i_m, i_o2, i_o3) & 
  * V2_(s_o4, s_o1, s_v1)%array(i_o4, i_o1, i_v1)
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

end subroutine g_sigma_ooov_ooov_no0_x36



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x36(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x36(sa, ia, xaaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x36



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x36(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1, s_m, i_m
! S2(i,k,m,a)  <-- 
! (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o3,o2,m,o1,o4,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(IEOR(s_o3,s_o2),s_m) == IEOR(IEOR(s_o1,s_o4),s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * X_(s_o4, s_o1, s_m, s_o2, s_o3)%array(i_o4, i_o1, i_m, i_o2, i_o3)
end do ! Irrep Loop
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
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no1_x36



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x37(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x37(sa, ia, sv1, iv1, av2_i, h2_i, xaaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x37



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x37(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_o4, i_o4, s_m, i_m, s_o3, i_o3
! X(o1,o2,o4,m,o3,a)  <-- 
! (    1.00000000)  T2(o1,o2,o4,v1) V2(a,v1,m,o3) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_m = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(IEOR(s_o1,s_o2),s_o4) == IEOR(IEOR(s_m,s_o3),s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_o4,s_v1) .and. &
IEOR(s_a,s_v1) == IEOR(s_m,s_o3)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_m, s_o4, s_o2, s_o1)%array(i_o3, i_m, i_o4, i_o2, i_o1) =  &
    X_(s_o3, s_m, s_o4, s_o2, s_o1)%array(i_o3, i_m, i_o4, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_o4, s_o2, s_o1)%array(i_o4, i_o2, i_o1) & 
  * V2_(s_o3, s_m, s_v1)%array(i_o3, i_m, i_v1)
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

end subroutine g_sigma_ooov_ooov_no0_x37



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x37(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x37(sa, ia, xaaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x37



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x37(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1, s_m, i_m
! S2(i,k,m,a)  <-- 
! (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o1,o2,o4,m,o3,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(IEOR(s_o1,s_o2),s_o4) == IEOR(IEOR(s_m,s_o3),s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * X_(s_o3, s_m, s_o4, s_o2, s_o1)%array(i_o3, i_m, i_o4, i_o2, i_o1)
end do ! Irrep Loop
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
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no1_x37



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x38(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x38(sa, ia, sv1, iv1, av2_i, h2_i, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x38



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x38(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o3, i_o3, s_o1, i_o1, s_o4, i_o4
! X(o2,o3,o4,a)  <-- 
! (    1.00000000)  T2(o2,o3,o1,v1) V2(a,v1,o1,o4) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o2,s_o3) == IEOR(s_o4,s_a) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_o1,s_v1) .and. &
IEOR(s_a,s_v1) == IEOR(s_o1,s_o4)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o3, s_o2)%array(i_o4, i_o3, i_o2) =  &
    X_(s_o4, s_o3, s_o2)%array(i_o4, i_o3, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2) & 
  * V2_(s_o4, s_o1, s_v1)%array(i_o4, i_o1, i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x38



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x38(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x38(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x38



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x38(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o3, i_o3, s_o4, i_o4
integer :: s_o2, i_o2
! S2(i,k,m,a)  <-- 
! (    1.00000000) D3(i,m,k,o3,o4,o2) X(o2,o3,o4,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o3,s_o4),s_o2) .and. &
IEOR(s_o2,s_o3) == IEOR(s_o4,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o4, s_o3, s_k, s_m, s_i)%array(i_o2, i_o4, i_o3, i_k, i_m, i_i) & 
  * X_(s_o4, s_o3, s_o2)%array(i_o4, i_o3, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x38



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x39(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x39(sa, ia, sv1, iv1, av2_i, h2_i, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x39



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x39(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_o1, i_o1, s_m, i_m
! X(o3,o2,m,a)  <-- 
! (    1.00000000)  T2(o3,o2,o1,v1) V2(a,v1,m,o1) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o3,s_o2) == IEOR(s_m,s_a) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_o1,s_v1) .and. &
IEOR(s_a,s_v1) == IEOR(s_m,s_o1)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o2, s_o3)%array(i_m, i_o2, i_o3) =  &
    X_(s_m, s_o2, s_o3)%array(i_m, i_o2, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o2, s_o3)%array(i_o1, i_o2, i_o3) & 
  * V2_(s_o1, s_m, s_v1)%array(i_o1, i_m, i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x39



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x39(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x39(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x39



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x39(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_m, i_m
! S2(i,k,m,a)  <-- 
! (    1.00000000) D2(i,o3,k,o2) X(o3,o2,m,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o3,s_o2) == IEOR(s_m,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * X_(s_m, s_o2, s_o3)%array(i_m, i_o2, i_o3)
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

end subroutine g_sigma_ooov_ooov_no1_x39



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x40(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x40(sa, ia, sv1, iv1, av2_i, h2_i, xaaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x40



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x40(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o4, i_o4, s_o5, i_o5, s_o1, i_o1, s_o3, i_o3, s_o2, i_o2
! X(o4,o5,o1,o3,o2,a)  <-- 
! (    1.00000000)  T2(o4,o5,o1,v1) V2(a,o3,o2,v1) 
do s_o4 = 0, nir-1
do s_o5 = 0, nir-1
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(IEOR(s_o4,s_o5),s_o1) == IEOR(IEOR(s_o3,s_o2),s_a) .and. & 
IEOR(s_o4,s_o5) == IEOR(s_o1,s_v1) .and. &
IEOR(s_a,s_o3) == IEOR(s_o2,s_v1)) then
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_(s_o2, s_o3, s_o1, s_o5, s_o4)%array(i_o2, i_o3, i_o1, i_o5, i_o4) =  &
    X_(s_o2, s_o3, s_o1, s_o5, s_o4)%array(i_o2, i_o3, i_o1, i_o5, i_o4) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o5, s_o4)%array(i_o1, i_o5, i_o4) & 
  * V2_(s_v1, s_o2, s_o3)%array(i_v1, i_o2, i_o3)
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

end subroutine g_sigma_ooov_ooov_no0_x40



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x40(sa, ia, si, ii, sm, im, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, si, ii, sm, im
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x40(sa, ia, si, ii, sm, im, xaaaaa, av2_i2, d4_ij, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x40



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x40(s_a, i_a, s_i, i_i, s_m, i_m, X_, S2_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_m, s_m
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_o3, i_o3, s_o1, i_o1, s_o4, i_o4, s_o2, i_o2
integer :: s_o5, i_o5
! S2(i,k,m,a)  <-- 
! (    1.00000000) D4(i,m,k,o3,o1,o4,o2,o5) X(o4,o5,o1,o3,o2,a) 
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
do s_o5 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),IEOR(s_k,s_o3)) == IEOR(IEOR(s_o1,s_o4),IEOR(s_o2,s_o5)) .and. &
IEOR(IEOR(s_o4,s_o5),s_o1) == IEOR(IEOR(s_o3,s_o2),s_a)) then
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D4_(s_o5, s_o2, s_o4, s_o1, s_o3, s_k)%array(i_o5, i_o2, i_o4, i_o1, i_o3, i_k) & 
  * X_(s_o2, s_o3, s_o1, s_o5, s_o4)%array(i_o2, i_o3, i_o1, i_o5, i_o4)
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

end subroutine g_sigma_ooov_ooov_no1_x40



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x41(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x41(sa, ia, sv1, iv1, av2_i, h2_i, xaaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x41



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x41(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o1, i_o1, s_m, i_m, s_o2, i_o2, s_o4, i_o4
! X(o3,o1,m,o2,o4,a)  <-- 
! (    1.00000000)  T2(o3,o1,m,v1) V2(a,o2,o4,v1) 
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(IEOR(s_o3,s_o1),s_m) == IEOR(IEOR(s_o2,s_o4),s_a) .and. & 
IEOR(s_o3,s_o1) == IEOR(s_m,s_v1) .and. &
IEOR(s_a,s_o2) == IEOR(s_o4,s_v1)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o2, s_m, s_o1, s_o3)%array(i_o4, i_o2, i_m, i_o1, i_o3) =  &
    X_(s_o4, s_o2, s_m, s_o1, s_o3)%array(i_o4, i_o2, i_m, i_o1, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_m, s_o1, s_o3)%array(i_m, i_o1, i_o3) & 
  * V2_(s_v1, s_o4, s_o2)%array(i_v1, i_o4, i_o2)
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

end subroutine g_sigma_ooov_ooov_no0_x41



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x41(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x41(sa, ia, xaaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x41



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x41(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1, s_m, i_m
! S2(i,k,m,a)  <-- 
! (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o3,o1,m,o2,o4,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(IEOR(s_o3,s_o1),s_m) == IEOR(IEOR(s_o2,s_o4),s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * X_(s_o4, s_o2, s_m, s_o1, s_o3)%array(i_o4, i_o2, i_m, i_o1, i_o3)
end do ! Irrep Loop
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
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no1_x41



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x42(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x42(sa, ia, sv1, iv1, av2_i, h2_i, xaaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x42



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x42(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3, s_o4, i_o4, s_o2, i_o2, s_m, i_m
! X(o1,o3,o4,o2,m,a)  <-- 
! (    1.00000000)  T2(o1,o3,o4,v1) V2(a,o2,m,v1) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(IEOR(s_o1,s_o3),s_o4) == IEOR(IEOR(s_o2,s_m),s_a) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o4,s_v1) .and. &
IEOR(s_a,s_o2) == IEOR(s_m,s_v1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o2, s_o4, s_o3, s_o1)%array(i_m, i_o2, i_o4, i_o3, i_o1) =  &
    X_(s_m, s_o2, s_o4, s_o3, s_o1)%array(i_m, i_o2, i_o4, i_o3, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_o4, s_o3, s_o1)%array(i_o4, i_o3, i_o1) & 
  * V2_(s_v1, s_m, s_o2)%array(i_v1, i_m, i_o2)
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

end subroutine g_sigma_ooov_ooov_no0_x42



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x42(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x42(sa, ia, xaaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x42



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x42(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1, s_m, i_m
! S2(i,k,m,a)  <-- 
! (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o1,o3,o4,o2,m,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(IEOR(s_o1,s_o3),s_o4) == IEOR(IEOR(s_o2,s_m),s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * X_(s_m, s_o2, s_o4, s_o3, s_o1)%array(i_m, i_o2, i_o4, i_o3, i_o1)
end do ! Irrep Loop
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
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no1_x42



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x43(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x43(sa, ia, sv1, iv1, av2_i, h2_i, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x43



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x43(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o3, i_o3, s_o4, i_o4, s_o1, i_o1
! X(o2,o3,o1,a)  <-- 
! (    1.00000000)  T2(o2,o3,o4,v1) V2(a,o4,o1,v1) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_o2,s_o3) == IEOR(s_o1,s_a) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_o4,s_v1) .and. &
IEOR(s_a,s_o4) == IEOR(s_o1,s_v1)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2) =  &
    X_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_o4, s_o3, s_o2)%array(i_o4, i_o3, i_o2) & 
  * V2_(s_v1, s_o1, s_o4)%array(i_v1, i_o1, i_o4)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x43



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x43(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x43(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x43



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x43(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o2, i_o2, s_o1, i_o1
integer :: s_o3, i_o3
! S2(i,k,m,a)  <-- 
! (    1.00000000) D3(i,m,k,o2,o1,o3) X(o2,o3,o1,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o2,s_o1),s_o3) .and. &
IEOR(s_o2,s_o3) == IEOR(s_o1,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o3, s_o1, s_o2, s_k, s_m, s_i)%array(i_o3, i_o1, i_o2, i_k, i_m, i_i) & 
  * X_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x43



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x44(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x44(sa, ia, sv1, iv1, av2_i, h2_i, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x44



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x44(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_o3, i_o3, s_m, i_m
! X(o2,o1,m,a)  <-- 
! (    1.00000000)  T2(o2,o1,o3,v1) V2(a,o3,m,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_m,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_o3,s_v1) .and. &
IEOR(s_a,s_o3) == IEOR(s_m,s_v1)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) =  &
    X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_o3, s_o1, s_o2)%array(i_o3, i_o1, i_o2) & 
  * V2_(s_v1, s_m, s_o3)%array(i_v1, i_m, i_o3)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ooov_no0_x44



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x44(sa, ia, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x44(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x44



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no1_x44(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2, s_m, i_m
! S2(i,k,m,a)  <-- 
! (    1.00000000) D2(i,o1,k,o2) X(o2,o1,m,a) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_o2) .and. &
IEOR(s_o2,s_o1) == IEOR(s_m,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i) & 
  * X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2)
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

end subroutine g_sigma_ooov_ooov_no1_x44



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x45(sa, ia, Ecas, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Ecas
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no0_x45(sa, ia, Ecas, av2_i, av2_i2, d3, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no0_x45



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x45(s_a, i_a, Ecas, T2_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Ecas
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o2, i_o2, s_o1, i_o1
integer :: s_o3, i_o3
! S2(i,k,m,a)  <-- 
! (    1.00000000) Ecas D3(i,m,k,o2,o1,o3) T2(o3,o2,o1,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o2,s_o1),s_o3) .and. &
IEOR(s_o3,s_o2) == IEOR(s_o1,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * Ecas & 
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

end subroutine g_sigma_ooov_ooov_no0_x45



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x46(sa, ia, Ecas, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Ecas
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no0_x46(sa, ia, Ecas, av2_i, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no0_x46



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ooov_no0_x46(s_a, i_a, Ecas, T2_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Ecas
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2, s_m, i_m
! S2(i,k,m,a)  <-- 
! (    1.00000000) Ecas D2(i,o1,k,o2) T2(o1,o2,m,a) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_o2) .and. &
IEOR(s_o1,s_o2) == IEOR(s_m,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * Ecas & 
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

end subroutine g_sigma_ooov_ooov_no0_x46

