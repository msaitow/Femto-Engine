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
subroutine g_if_sigma_ooov_ccvv_no0_x0(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

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

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xa(sleft, x, nir, nsym, psym) ! -> xa (allocate) 
call g_sigma_ooov_ccvv_no0_x0(sa, ia, sv1, iv1, av2_i, h2_i, xa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xa)

end subroutine g_if_sigma_ooov_ccvv_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ccvv_no0_x0(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_c1, i_c1, s_o1, i_o1
! X(o1,a)  <-- 
! (    1.00000000)  T2(c2,c1,v1,a) V2(v1,c2,c1,o1) 
do s_c2 = 0, nir-1
do s_c1 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_o1,s_a) == 0 .and. & 
IEOR(s_c2,s_c1) == IEOR(s_v1,s_a) .and. &
IEOR(s_v1,s_c2) == IEOR(s_c1,s_o1)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1)%array(i_o1) =  &
    X_(s_o1)%array(i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_c1, s_c2)%array(i_v1, i_c1, i_c2) & 
  * V2_(s_o1, s_c1, s_c2)%array(i_o1, i_c1, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ccvv_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ccvv_no1_x0(sa, ia, X, S2, nir, nsym, psym)

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

call set_symblock_xa(sleft, x, nir, nsym, psym) ! -> xa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ccvv_no1_x0(sa, ia, xa, av2_i2, d2, nir, nsym, psym)

deallocate(xa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ccvv_no1_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ccvv_no1_x0(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! S2(i,k,m,a)  <-- 
! (   -2.00000000) D2(i,m,k,o1) X(o1,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_m) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 2.00000000d+00 & 
  * D2_(s_o1, s_k, s_m, s_i)%array(i_o1, i_k, i_m, i_i) & 
  * X_(s_o1)%array(i_o1)
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

end subroutine g_sigma_ooov_ccvv_no1_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ccvv_no0_x1(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

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

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xa(sleft, x, nir, nsym, psym) ! -> xa (allocate) 
call g_sigma_ooov_ccvv_no0_x1(sa, ia, sv1, iv1, av2_i, h2_i, xa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xa)

end subroutine g_if_sigma_ooov_ccvv_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ccvv_no0_x1(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_c1, i_c1, s_o1, i_o1
! X(o1,a)  <-- 
! (    1.00000000)  T2(c2,c1,v1,a) V2(v1,c1,c2,o1) 
do s_c2 = 0, nir-1
do s_c1 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_o1,s_a) == 0 .and. & 
IEOR(s_c2,s_c1) == IEOR(s_v1,s_a) .and. &
IEOR(s_v1,s_c1) == IEOR(s_c2,s_o1)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1)%array(i_o1) =  &
    X_(s_o1)%array(i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_c1, s_c2)%array(i_v1, i_c1, i_c2) & 
  * V2_(s_o1, s_c2, s_c1)%array(i_o1, i_c2, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ccvv_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ccvv_no1_x1(sa, ia, X, S2, nir, nsym, psym)

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

call set_symblock_xa(sleft, x, nir, nsym, psym) ! -> xa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ccvv_no1_x1(sa, ia, xa, av2_i2, d2, nir, nsym, psym)

deallocate(xa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ccvv_no1_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ccvv_no1_x1(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! S2(i,k,m,a)  <-- 
! (    1.00000000) D2(i,m,k,o1) X(o1,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_m) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_m, s_i)%array(i_o1, i_k, i_m, i_i) & 
  * X_(s_o1)%array(i_o1)
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

end subroutine g_sigma_ooov_ccvv_no1_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ccvv_no0_x2(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

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
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xa(sleft, x, nir, nsym, psym) ! -> xa (allocate) 
call g_sigma_ooov_ccvv_no0_x2(sa, ia, sv1, iv1, av2_i, h2_i, xa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xa)

end subroutine g_if_sigma_ooov_ccvv_no0_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ccvv_no0_x2(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_c1, i_c1, s_o1, i_o1
! X(o1,a)  <-- 
! (    1.00000000)  T2(c2,c1,a,v1) V2(v1,c2,c1,o1) 
do s_c2 = 0, nir-1
do s_c1 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_o1,s_a) == 0 .and. & 
IEOR(s_c2,s_c1) == IEOR(s_a,s_v1) .and. &
IEOR(s_v1,s_c2) == IEOR(s_c1,s_o1)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1)%array(i_o1) =  &
    X_(s_o1)%array(i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_c1, s_c2)%array(i_a, i_c1, i_c2) & 
  * V2_(s_o1, s_c1, s_c2)%array(i_o1, i_c1, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ccvv_no0_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ccvv_no1_x2(sa, ia, X, S2, nir, nsym, psym)

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

call set_symblock_xa(sleft, x, nir, nsym, psym) ! -> xa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ccvv_no1_x2(sa, ia, xa, av2_i2, d2, nir, nsym, psym)

deallocate(xa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ccvv_no1_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ccvv_no1_x2(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! S2(i,k,m,a)  <-- 
! (    1.00000000) D2(i,m,k,o1) X(o1,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_m) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_m, s_i)%array(i_o1, i_k, i_m, i_i) & 
  * X_(s_o1)%array(i_o1)
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

end subroutine g_sigma_ooov_ccvv_no1_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ccvv_no0_x3(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

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
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xa(sleft, x, nir, nsym, psym) ! -> xa (allocate) 
call g_sigma_ooov_ccvv_no0_x3(sa, ia, sv1, iv1, av2_i, h2_i, xa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xa)

end subroutine g_if_sigma_ooov_ccvv_no0_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ccvv_no0_x3(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_c2, i_c2, s_o1, i_o1
! X(o1,a)  <-- 
! (    1.00000000)  T2(c1,c2,a,v1) V2(v1,c2,c1,o1) 
do s_c1 = 0, nir-1
do s_c2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_o1,s_a) == 0 .and. & 
IEOR(s_c1,s_c2) == IEOR(s_a,s_v1) .and. &
IEOR(s_v1,s_c2) == IEOR(s_c1,s_o1)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1)%array(i_o1) =  &
    X_(s_o1)%array(i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_c2, s_c1)%array(i_a, i_c2, i_c1) & 
  * V2_(s_o1, s_c1, s_c2)%array(i_o1, i_c1, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ccvv_no0_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ccvv_no1_x3(sa, ia, X, S2, nir, nsym, psym)

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

call set_symblock_xa(sleft, x, nir, nsym, psym) ! -> xa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ccvv_no1_x3(sa, ia, xa, av2_i2, d2, nir, nsym, psym)

deallocate(xa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ccvv_no1_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_ccvv_no1_x3(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! S2(i,k,m,a)  <-- 
! (   -2.00000000) D2(i,m,k,o1) X(o1,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_m) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 2.00000000d+00 & 
  * D2_(s_o1, s_k, s_m, s_i)%array(i_o1, i_k, i_m, i_i) & 
  * X_(s_o1)%array(i_o1)
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

end subroutine g_sigma_ooov_ccvv_no1_x3

