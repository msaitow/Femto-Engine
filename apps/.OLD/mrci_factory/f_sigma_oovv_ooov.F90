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
subroutine g_if_sigma_oovv_ooov_no0_x0(sa, ia, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_oovv_ooov_no0_x0(sa, ia, av2_i, xaaa, d3, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x0(s_a, i_a, T2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! X(i,k,o2,a)  <-- 
! (    1.00000000)  D3(i,o3,k,o2,o4,o1) T2(o1,o3,o4,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_o2,s_a) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o4,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o2, s_k, s_i)%array(i_o2, i_k, i_i) =  &
    X_(s_o2, s_k, s_i)%array(i_o2, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * T2_(s_o4, s_o3, s_o1)%array(i_o4, i_o3, i_o1)
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

end subroutine g_sigma_oovv_ooov_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x0(sa, ia, sc, ic, X, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x0(sa, ia, sc, ic, xaaa, h1, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x0(s_a, i_a, s_c, i_c, X_, h_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o2, i_o2
! S2(i,k,a,c)  <-- 
! (    1.00000000) X(i,k,o2,a) h(o2,c) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_k) == IEOR(s_o2,s_a) .and. &
IEOR(s_o2,s_c) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * X_(s_o2, s_k, s_i)%array(i_o2, i_k, i_i) & 
  * h_(s_c, s_o2)%array(i_c, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no1_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x1(sc, ic, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_oovv_ooov_no0_x1(sc, ic, av2_i, xaaa, d3, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x1(s_c, i_c, T2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o3, i_o3, s_o1, i_o1
integer :: s_o4, i_o4
! X(i,o2,k,c)  <-- 
! (    1.00000000)  D3(i,o2,k,o3,o1,o4) T2(o4,o3,o1,c) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_i,s_o2) == IEOR(s_k,s_c) .and. & 
IEOR(IEOR(s_i,s_o2),s_k) == IEOR(IEOR(s_o3,s_o1),s_o4) .and. &
IEOR(s_o4,s_o3) == IEOR(s_o1,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_k, s_o2, s_i)%array(i_k, i_o2, i_i) =  &
    X_(s_k, s_o2, s_i)%array(i_k, i_o2, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o4, s_o1, s_o3, s_k, s_o2, s_i)%array(i_o4, i_o1, i_o3, i_k, i_o2, i_i) & 
  * T2_(s_o1, s_o3, s_o4)%array(i_o1, i_o3, i_o4)
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

end subroutine g_sigma_oovv_ooov_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x1(sc, ic, X, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sc

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x1(sc, ic, xaaa, h1, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x1(s_c, i_c, X_, h_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_a, i_a
! S2(i,k,a,c)  <-- 
! (    1.00000000) X(i,o2,k,c) h(o2,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_c) .and. &
IEOR(s_o2,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * X_(s_k, s_o2, s_i)%array(i_k, i_o2, i_i) & 
  * h_(s_a, s_o2)%array(i_a, i_o2)
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

end subroutine g_sigma_oovv_ooov_no1_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x2(sa, ia, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_oovv_ooov_no0_x2(sa, ia, av2_i, xaaa, d2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x2(s_a, i_a, T2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o1, i_o1
! X(i,k,o1,a)  <-- 
! (    1.00000000)  D2(i,o3,k,o2) T2(o2,o3,o1,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_o1,s_a) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o2,s_o3) == IEOR(s_o1,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) =  &
    X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * T2_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2)
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

end subroutine g_sigma_oovv_ooov_no0_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x2(sa, ia, sc, ic, X, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x2(sa, ia, sc, ic, xaaa, h1, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x2(s_a, i_a, s_c, i_c, X_, h_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1
! S2(i,k,a,c)  <-- 
! (    1.00000000) X(i,k,o1,a) h(o1,c) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_k) == IEOR(s_o1,s_a) .and. &
IEOR(s_o1,s_c) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) & 
  * h_(s_c, s_o1)%array(i_c, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no1_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x3(sc, ic, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_oovv_ooov_no0_x3(sc, ic, av2_i, xaaa, d2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x3(s_c, i_c, T2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2, s_o3, i_o3
! X(i,k,o3,c)  <-- 
! (    1.00000000)  D2(i,o1,k,o2) T2(o1,o2,o3,c) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_o3,s_c) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_o2) .and. &
IEOR(s_o1,s_o2) == IEOR(s_o3,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_k, s_i)%array(i_o3, i_k, i_i) =  &
    X_(s_o3, s_k, s_i)%array(i_o3, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i) & 
  * T2_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1)
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

end subroutine g_sigma_oovv_ooov_no0_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x3(sc, ic, X, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sc

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x3(sc, ic, xaaa, h1, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x3(s_c, i_c, X_, h_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o3, i_o3, s_a, i_a
! S2(i,k,a,c)  <-- 
! (    1.00000000) X(i,k,o3,c) h(o3,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_k) == IEOR(s_o3,s_c) .and. &
IEOR(s_o3,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * X_(s_o3, s_k, s_i)%array(i_o3, i_k, i_i) & 
  * h_(s_a, s_o3)%array(i_a, i_o3)
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

end subroutine g_sigma_oovv_ooov_no1_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_y4(sc1, ic1, V2, Y0, nir, nsym, psym)

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
call set_symblock_Yav(sleft2, Y0, nir, nsym, psym) ! -> Yav (allocate) 
call g_sigma_oovv_ooov_y4(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_oovv_ooov_y4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_y4(s_c1, i_c1, V2_, Y0_, nir, nsym, psym)

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

integer :: s_o2, i_o2, s_c, i_c
! Y0(o2,c)  <-- 
! (    1.00000000)  V2(c1,c1,o2,c) 
do s_o2 = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_o2,s_c) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o2,s_c)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y0_(s_c, s_o2)%array(i_c, i_o2) =  &
    Y0_(s_c, s_o2)%array(i_c, i_o2) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_o2, s_c1)%array(i_c, i_o2, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_y4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x4(sa, ia, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_oovv_ooov_no0_x4(sa, ia, av2_i, xaaa, d3, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x4(s_a, i_a, T2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! X(i,k,o2,a)  <-- 
! (    1.00000000)  D3(i,o3,k,o2,o4,o1) T2(o1,o3,o4,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_o2,s_a) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o4,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o2, s_k, s_i)%array(i_o2, i_k, i_i) =  &
    X_(s_o2, s_k, s_i)%array(i_o2, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * T2_(s_o4, s_o3, s_o1)%array(i_o4, i_o3, i_o1)
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

end subroutine g_sigma_oovv_ooov_no0_x4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x4(sa, ia, sc, ic, X, Y0, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X(*), Y0(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
sleft2 = 0
call set_symblock_Yav(sleft2, Y0, nir, nsym, psym) ! -> Yav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x4(sa, ia, sc, ic, xaaa, Yav, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(Yav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x4(s_a, i_a, s_c, i_c, X_, Y0_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y0_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o2, i_o2
! S2(i,k,a,c)  <-- 
! (    2.00000000) X(i,k,o2,a) Y0(o2,c) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_k) == IEOR(s_o2,s_a) .and. &
IEOR(s_o2,s_c) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * X_(s_o2, s_k, s_i)%array(i_o2, i_k, i_i) & 
  * Y0_(s_c, s_o2)%array(i_c, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no1_x4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_y5(sc1, ic1, V2, Y1, nir, nsym, psym)

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
call set_symblock_Yav(sleft2, Y1, nir, nsym, psym) ! -> Yav (allocate) 
call g_sigma_oovv_ooov_y5(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_oovv_ooov_y5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_y5(s_c1, i_c1, V2_, Y1_, nir, nsym, psym)

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

integer :: s_o2, i_o2, s_c, i_c
! Y1(o2,c)  <-- 
! (    1.00000000)  V2(c1,o2,c1,c) 
do s_o2 = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_o2,s_c) == 0 .and. & 
IEOR(s_c1,s_o2) == IEOR(s_c1,s_c)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y1_(s_c, s_o2)%array(i_c, i_o2) =  &
    Y1_(s_c, s_o2)%array(i_c, i_o2) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_c1, s_o2)%array(i_c, i_c1, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_y5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x5(sa, ia, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_oovv_ooov_no0_x5(sa, ia, av2_i, xaaa, d3, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x5(s_a, i_a, T2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! X(i,k,o2,a)  <-- 
! (    1.00000000)  D3(i,o3,k,o2,o4,o1) T2(o1,o3,o4,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_o2,s_a) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o4,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o2, s_k, s_i)%array(i_o2, i_k, i_i) =  &
    X_(s_o2, s_k, s_i)%array(i_o2, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * T2_(s_o4, s_o3, s_o1)%array(i_o4, i_o3, i_o1)
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

end subroutine g_sigma_oovv_ooov_no0_x5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x5(sa, ia, sc, ic, X, Y1, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X(*), Y1(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
sleft2 = 0
call set_symblock_Yav(sleft2, Y1, nir, nsym, psym) ! -> Yav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x5(sa, ia, sc, ic, xaaa, Yav, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(Yav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x5(s_a, i_a, s_c, i_c, X_, Y1_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y1_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o2, i_o2
! S2(i,k,a,c)  <-- 
! (   -1.00000000) X(i,k,o2,a) Y1(o2,c) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_k) == IEOR(s_o2,s_a) .and. &
IEOR(s_o2,s_c) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * X_(s_o2, s_k, s_i)%array(i_o2, i_k, i_i) & 
  * Y1_(s_c, s_o2)%array(i_c, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no1_x5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_y6(sc1, ic1, V2, Y2, nir, nsym, psym)

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
call set_symblock_Yav(sleft2, Y2, nir, nsym, psym) ! -> Yav (allocate) 
call g_sigma_oovv_ooov_y6(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_oovv_ooov_y6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_y6(s_c1, i_c1, V2_, Y2_, nir, nsym, psym)

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

integer :: s_o3, i_o3, s_a, i_a
! Y2(o3,a)  <-- 
! (    1.00000000)  V2(c1,c1,o3,a) 
do s_o3 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o3,s_a) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o3,s_a)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y2_(s_a, s_o3)%array(i_a, i_o3) =  &
    Y2_(s_a, s_o3)%array(i_a, i_o3) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_o3, s_c1)%array(i_a, i_o3, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_y6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x6(sc, ic, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_oovv_ooov_no0_x6(sc, ic, av2_i, xaaa, d3, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x6(s_c, i_c, T2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! X(i,o3,k,c)  <-- 
! (    1.00000000)  D3(i,o3,k,o2,o4,o1) T2(o1,o2,o4,c) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_o3) == IEOR(s_k,s_c) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_o4,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_k, s_o3, s_i)%array(i_k, i_o3, i_i) =  &
    X_(s_k, s_o3, s_i)%array(i_k, i_o3, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * T2_(s_o4, s_o2, s_o1)%array(i_o4, i_o2, i_o1)
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

end subroutine g_sigma_oovv_ooov_no0_x6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x6(sc, ic, X, Y2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X(*), Y2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sc

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
sleft2 = 0
call set_symblock_Yav(sleft2, Y2, nir, nsym, psym) ! -> Yav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x6(sc, ic, xaaa, Yav, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(Yav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x6(s_c, i_c, X_, Y2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y2_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_a, i_a
! S2(i,k,a,c)  <-- 
! (    2.00000000) X(i,o3,k,c) Y2(o3,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_c) .and. &
IEOR(s_o3,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * X_(s_k, s_o3, s_i)%array(i_k, i_o3, i_i) & 
  * Y2_(s_a, s_o3)%array(i_a, i_o3)
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

end subroutine g_sigma_oovv_ooov_no1_x6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_y7(sc1, ic1, V2, Y3, nir, nsym, psym)

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
call set_symblock_Yav(sleft2, Y3, nir, nsym, psym) ! -> Yav (allocate) 
call g_sigma_oovv_ooov_y7(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_oovv_ooov_y7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_y7(s_c1, i_c1, V2_, Y3_, nir, nsym, psym)

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

integer :: s_o3, i_o3, s_a, i_a
! Y3(o3,a)  <-- 
! (    1.00000000)  V2(c1,o3,c1,a) 
do s_o3 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o3,s_a) == 0 .and. & 
IEOR(s_c1,s_o3) == IEOR(s_c1,s_a)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y3_(s_a, s_o3)%array(i_a, i_o3) =  &
    Y3_(s_a, s_o3)%array(i_a, i_o3) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_c1, s_o3)%array(i_a, i_c1, i_o3)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_y7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x7(sc, ic, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_oovv_ooov_no0_x7(sc, ic, av2_i, xaaa, d3, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x7(s_c, i_c, T2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! X(i,o3,k,c)  <-- 
! (    1.00000000)  D3(i,o3,k,o2,o4,o1) T2(o1,o2,o4,c) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_o3) == IEOR(s_k,s_c) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_o4,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_k, s_o3, s_i)%array(i_k, i_o3, i_i) =  &
    X_(s_k, s_o3, s_i)%array(i_k, i_o3, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * T2_(s_o4, s_o2, s_o1)%array(i_o4, i_o2, i_o1)
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

end subroutine g_sigma_oovv_ooov_no0_x7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x7(sc, ic, X, Y3, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X(*), Y3(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sc

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
sleft2 = 0
call set_symblock_Yav(sleft2, Y3, nir, nsym, psym) ! -> Yav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x7(sc, ic, xaaa, Yav, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(Yav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x7(s_c, i_c, X_, Y3_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y3_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_a, i_a
! S2(i,k,a,c)  <-- 
! (   -1.00000000) X(i,o3,k,c) Y3(o3,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_c) .and. &
IEOR(s_o3,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * X_(s_k, s_o3, s_i)%array(i_k, i_o3, i_i) & 
  * Y3_(s_a, s_o3)%array(i_a, i_o3)
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

end subroutine g_sigma_oovv_ooov_no1_x7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_y8(sc1, ic1, V2, Y4, nir, nsym, psym)

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
call set_symblock_Yav(sleft2, Y4, nir, nsym, psym) ! -> Yav (allocate) 
call g_sigma_oovv_ooov_y8(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_oovv_ooov_y8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_y8(s_c1, i_c1, V2_, Y4_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_c, i_c
! Y4(o1,c)  <-- 
! (    1.00000000)  V2(c1,c1,o1,c) 
do s_o1 = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_o1,s_c) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_c)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y4_(s_c, s_o1)%array(i_c, i_o1) =  &
    Y4_(s_c, s_o1)%array(i_c, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_o1, s_c1)%array(i_c, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_y8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x8(sa, ia, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_oovv_ooov_no0_x8(sa, ia, av2_i, xaaa, d2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x8(s_a, i_a, T2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o1, i_o1
! X(i,k,o1,a)  <-- 
! (    1.00000000)  D2(i,o3,k,o2) T2(o2,o3,o1,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_o1,s_a) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o2,s_o3) == IEOR(s_o1,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) =  &
    X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * T2_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2)
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

end subroutine g_sigma_oovv_ooov_no0_x8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x8(sa, ia, sc, ic, X, Y4, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X(*), Y4(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
sleft2 = 0
call set_symblock_Yav(sleft2, Y4, nir, nsym, psym) ! -> Yav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x8(sa, ia, sc, ic, xaaa, Yav, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(Yav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x8(s_a, i_a, s_c, i_c, X_, Y4_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y4_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1
! S2(i,k,a,c)  <-- 
! (    2.00000000) X(i,k,o1,a) Y4(o1,c) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_k) == IEOR(s_o1,s_a) .and. &
IEOR(s_o1,s_c) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) & 
  * Y4_(s_c, s_o1)%array(i_c, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no1_x8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_y9(sc1, ic1, V2, Y5, nir, nsym, psym)

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
call set_symblock_Yav(sleft2, Y5, nir, nsym, psym) ! -> Yav (allocate) 
call g_sigma_oovv_ooov_y9(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_oovv_ooov_y9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_y9(s_c1, i_c1, V2_, Y5_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_c, i_c
! Y5(o1,c)  <-- 
! (    1.00000000)  V2(c1,o1,c1,c) 
do s_o1 = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_o1,s_c) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_c)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y5_(s_c, s_o1)%array(i_c, i_o1) =  &
    Y5_(s_c, s_o1)%array(i_c, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_c1, s_o1)%array(i_c, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_y9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x9(sa, ia, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_oovv_ooov_no0_x9(sa, ia, av2_i, xaaa, d2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x9(s_a, i_a, T2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o1, i_o1
! X(i,k,o1,a)  <-- 
! (    1.00000000)  D2(i,o3,k,o2) T2(o2,o3,o1,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_o1,s_a) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o2,s_o3) == IEOR(s_o1,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) =  &
    X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * T2_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2)
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

end subroutine g_sigma_oovv_ooov_no0_x9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x9(sa, ia, sc, ic, X, Y5, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X(*), Y5(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
sleft2 = 0
call set_symblock_Yav(sleft2, Y5, nir, nsym, psym) ! -> Yav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x9(sa, ia, sc, ic, xaaa, Yav, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(Yav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x9(s_a, i_a, s_c, i_c, X_, Y5_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y5_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1
! S2(i,k,a,c)  <-- 
! (   -1.00000000) X(i,k,o1,a) Y5(o1,c) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_k) == IEOR(s_o1,s_a) .and. &
IEOR(s_o1,s_c) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) & 
  * Y5_(s_c, s_o1)%array(i_c, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no1_x9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_y10(sc1, ic1, V2, Y6, nir, nsym, psym)

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
call set_symblock_Yav(sleft2, Y6, nir, nsym, psym) ! -> Yav (allocate) 
call g_sigma_oovv_ooov_y10(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_oovv_ooov_y10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_y10(s_c1, i_c1, V2_, Y6_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_a, i_a
! Y6(o1,a)  <-- 
! (    1.00000000)  V2(c1,c1,o1,a) 
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o1,s_a) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_a)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y6_(s_a, s_o1)%array(i_a, i_o1) =  &
    Y6_(s_a, s_o1)%array(i_a, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_o1, s_c1)%array(i_a, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_y10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x10(sc, ic, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_oovv_ooov_no0_x10(sc, ic, av2_i, xaaa, d2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x10(s_c, i_c, T2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o1, i_o1
! X(i,k,o1,c)  <-- 
! (    1.00000000)  D2(i,o3,k,o2) T2(o3,o2,o1,c) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_o1,s_c) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o3,s_o2) == IEOR(s_o1,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) =  &
    X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * T2_(s_o1, s_o2, s_o3)%array(i_o1, i_o2, i_o3)
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

end subroutine g_sigma_oovv_ooov_no0_x10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x10(sc, ic, X, Y6, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X(*), Y6(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sc

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
sleft2 = 0
call set_symblock_Yav(sleft2, Y6, nir, nsym, psym) ! -> Yav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x10(sc, ic, xaaa, Yav, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(Yav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x10(s_c, i_c, X_, Y6_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y6_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_a, i_a
! S2(i,k,a,c)  <-- 
! (    2.00000000) X(i,k,o1,c) Y6(o1,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_k) == IEOR(s_o1,s_c) .and. &
IEOR(s_o1,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) & 
  * Y6_(s_a, s_o1)%array(i_a, i_o1)
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

end subroutine g_sigma_oovv_ooov_no1_x10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_y11(sc1, ic1, V2, Y7, nir, nsym, psym)

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
call set_symblock_Yav(sleft2, Y7, nir, nsym, psym) ! -> Yav (allocate) 
call g_sigma_oovv_ooov_y11(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_oovv_ooov_y11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_y11(s_c1, i_c1, V2_, Y7_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_a, i_a
! Y7(o1,a)  <-- 
! (    1.00000000)  V2(c1,o1,c1,a) 
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o1,s_a) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_a)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y7_(s_a, s_o1)%array(i_a, i_o1) =  &
    Y7_(s_a, s_o1)%array(i_a, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_c1, s_o1)%array(i_a, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_y11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x11(sc, ic, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_oovv_ooov_no0_x11(sc, ic, av2_i, xaaa, d2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x11(s_c, i_c, T2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o1, i_o1
! X(i,k,o1,c)  <-- 
! (    1.00000000)  D2(i,o3,k,o2) T2(o3,o2,o1,c) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_o1,s_c) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o3,s_o2) == IEOR(s_o1,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) =  &
    X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * T2_(s_o1, s_o2, s_o3)%array(i_o1, i_o2, i_o3)
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

end subroutine g_sigma_oovv_ooov_no0_x11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x11(sc, ic, X, Y7, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X(*), Y7(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sc

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
sleft2 = 0
call set_symblock_Yav(sleft2, Y7, nir, nsym, psym) ! -> Yav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x11(sc, ic, xaaa, Yav, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(Yav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x11(s_c, i_c, X_, Y7_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y7_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_a, i_a
! S2(i,k,a,c)  <-- 
! (   -1.00000000) X(i,k,o1,c) Y7(o1,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_k) == IEOR(s_o1,s_c) .and. &
IEOR(s_o1,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * X_(s_o1, s_k, s_i)%array(i_o1, i_k, i_i) & 
  * Y7_(s_a, s_o1)%array(i_a, i_o1)
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

end subroutine g_sigma_oovv_ooov_no1_x11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x12(sc, ic, si, ii, so4, io4, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, si, ii, so4, io4
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(si,IEOR(so4,sc))

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_oovv_ooov_no0_x12(sc, ic, si, ii, so4, io4, h2_i, xaaa, d4_ij, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x12(s_c, i_c, s_i, i_i, s_o4, i_o4, V2_, X_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_o4, s_o4
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_o3, i_o3, s_o5, i_o5, s_o1, i_o1, s_o6, i_o6
integer :: s_o2, i_o2
! X(i,o4,k,o6,o2,c)  <-- 
! (    1.00000000)  D4(i,o4,k,o3,o5,o1,o6,o2) V2(c,o3,o1,o5) 
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o5 = 0, nir-1
do s_o1 = 0, nir-1
do s_o6 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o6,s_o2),s_c) .and. & 
IEOR(IEOR(s_i,s_o4),IEOR(s_k,s_o3)) == IEOR(IEOR(s_o5,s_o1),IEOR(s_o6,s_o2)) .and. &
IEOR(s_c,s_o3) == IEOR(s_o1,s_o5)) then
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o6 = psym(I_BEGIN, I_O, s_o6), psym(I_END, I_O, s_o6)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_(s_o2, s_o6, s_k)%array(i_o2, i_o6, i_k) =  &
    X_(s_o2, s_o6, s_k)%array(i_o2, i_o6, i_k) &
  + 1.00000000d+00 & 
  * D4_(s_o2, s_o6, s_o1, s_o5, s_o3, s_k)%array(i_o2, i_o6, i_o1, i_o5, i_o3, i_k) & 
  * V2_(s_o5, s_o1, s_o3)%array(i_o5, i_o1, i_o3)
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

end subroutine g_sigma_oovv_ooov_no0_x12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x12(sa, ia, sc, ic, si, ii, so4, io4, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic, si, ii, so4, io4
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(si,IEOR(so4,sc))

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x12(sa, ia, sc, ic, si, ii, so4, io4, av2_i, xaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x12(s_a, i_a, s_c, i_c, s_i, i_i, s_o4, i_o4, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o4, s_o4
integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o6, i_o6, s_k, i_k
! S2(i,k,a,c)  <-- 
! (    1.00000000) T2(o2,o4,o6,a) X(i,o4,k,o6,o2,c) 
do s_o2 = 0, nir-1
do s_o6 = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o2,s_o4) == IEOR(s_o6,s_a) .and. &
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o6,s_o2),s_c)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o6 = psym(I_BEGIN, I_O, s_o6), psym(I_END, I_O, s_o6)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * T2_(s_o6, s_o4, s_o2)%array(i_o6, i_o4, i_o2) & 
  * X_(s_o2, s_o6, s_k)%array(i_o2, i_o6, i_k)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no1_x12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x13(sc, ic, so2, io2, so6, io6, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, so2, io2, so6, io6
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(so2, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(so2,IEOR(so6,sc))

call set_symblock_xaaaav(sleft, x, nir, nsym, psym) ! -> xaaaav (allocate) 
call g_sigma_oovv_ooov_no0_x13(sc, ic, so2, io2, so6, io6, av2_i, h2_i, xaaaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaav)

end subroutine g_if_sigma_oovv_ooov_no0_x13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x13(s_c, i_c, s_o2, i_o2, s_o6, i_o6, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_o2, s_o2
integer, intent(in) :: i_o6, s_o6
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o5, i_o5, s_o4, i_o4, s_o1, i_o1, s_o3, i_o3, s_a, i_a
! X(o5,o4,o1,o2,o6,o3,c,a)  <-- 
! (    1.00000000)  T2(o5,o4,o1,c) V2(o2,o6,o3,a) 
do s_o5 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(IEOR(s_o5,s_o4),IEOR(s_o1,s_o2)) == IEOR(IEOR(s_o6,s_o3),IEOR(s_c,s_a)) .and. & 
IEOR(s_o5,s_o4) == IEOR(s_o1,s_c) .and. &
IEOR(s_o2,s_o6) == IEOR(s_o3,s_a)) then
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o3, s_o1, s_o4, s_o5)%array(i_a, i_o3, i_o1, i_o4, i_o5) =  &
    X_(s_a, s_o3, s_o1, s_o4, s_o5)%array(i_a, i_o3, i_o1, i_o4, i_o5) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o4, s_o5)%array(i_o1, i_o4, i_o5) & 
  * V2_(s_a, s_o3, s_o6)%array(i_a, i_o3, i_o6)
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

end subroutine g_sigma_oovv_ooov_no0_x13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x13(sc, ic, so2, io2, so6, io6, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, so2, io2, so6, io6
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = IEOR(so2,IEOR(so6,sc))

call set_symblock_xaaaav(sleft, x, nir, nsym, psym) ! -> xaaaav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x13(sc, ic, so2, io2, so6, io6, xaaaav, av2_i2, d4_ij, nir, nsym, psym)

deallocate(xaaaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x13(s_c, i_c, s_o2, i_o2, s_o6, i_o6, X_, S2_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o2, s_o2
integer, intent(in) :: i_o6, s_o6
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o4, i_o4, s_o1, i_o1
integer :: s_o5, i_o5, s_a, i_a
! S2(i,k,a,c)  <-- 
! (    1.00000000) D4(o2,o6,i,o3,k,o4,o1,o5) X(o5,o4,o1,o2,o6,o3,c,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
do s_o5 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_o2,s_o6),IEOR(s_i,s_o3)) == IEOR(IEOR(s_k,s_o4),IEOR(s_o1,s_o5)) .and. &
IEOR(IEOR(s_o5,s_o4),IEOR(s_o1,s_o2)) == IEOR(IEOR(s_o6,s_o3),IEOR(s_c,s_a))) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D4_(s_o5, s_o1, s_o4, s_k, s_o3, s_i)%array(i_o5, i_o1, i_o4, i_k, i_o3, i_i) & 
  * X_(s_a, s_o3, s_o1, s_o4, s_o5)%array(i_a, i_o3, i_o1, i_o4, i_o5)
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

end subroutine g_sigma_oovv_ooov_no1_x13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x14(sc, ic, V2, X, nir, nsym, psym)

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

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_oovv_ooov_no0_x14(sc, ic, h2_i, xaaaaa, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x14(s_c, i_c, V2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o4, i_o4, s_k, i_k, s_o3, i_o3, s_o5, i_o5
integer :: s_o2, i_o2, s_o1, i_o1
! X(i,o4,k,o2,o1,c)  <-- 
! (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(c,o3,o1,o5) 
do s_i = 0, nir-1
do s_o4 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o5 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o2,s_o1),s_c) .and. & 
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_o5),s_o2) .and. &
IEOR(s_c,s_o3) == IEOR(s_o1,s_o5)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_o2, s_k, s_o4, s_i)%array(i_o1, i_o2, i_k, i_o4, i_i) =  &
    X_(s_o1, s_o2, s_k, s_o4, s_i)%array(i_o1, i_o2, i_k, i_o4, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o5, s_o3, s_k, s_o4, s_i)%array(i_o2, i_o5, i_o3, i_k, i_o4, i_i) & 
  * V2_(s_o5, s_o1, s_o3)%array(i_o5, i_o1, i_o3)
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

end subroutine g_sigma_oovv_ooov_no0_x14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x14(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x14(sa, ia, sc, ic, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x14(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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

integer :: s_o2, i_o2, s_o4, i_o4, s_o1, i_o1, s_i, i_i, s_k, i_k
! S2(i,k,a,c)  <-- 
! (    1.00000000) T2(o2,o4,o1,a) X(i,o4,k,o2,o1,c) 
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o2,s_o4) == IEOR(s_o1,s_a) .and. &
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o2,s_o1),s_c)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o4, s_o2)%array(i_o1, i_o4, i_o2) & 
  * X_(s_o1, s_o2, s_k, s_o4, s_i)%array(i_o1, i_o2, i_k, i_o4, i_i)
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

end subroutine g_sigma_oovv_ooov_no1_x14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x15(sa, ia, V2, X, nir, nsym, psym)

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

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_oovv_ooov_no0_x15(sa, ia, h2_i, xaaaaa, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x15(s_a, i_a, V2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o4, i_o4, s_k, i_k, s_o3, i_o3, s_o5, i_o5
integer :: s_o2, i_o2, s_o1, i_o1
! X(i,k,o3,o2,o1,a)  <-- 
! (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(a,o4,o1,o5) 
do s_i = 0, nir-1
do s_o4 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o5 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(IEOR(s_i,s_k),s_o3) == IEOR(IEOR(s_o2,s_o1),s_a) .and. & 
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_o5),s_o2) .and. &
IEOR(s_a,s_o4) == IEOR(s_o1,s_o5)) then
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

end subroutine g_sigma_oovv_ooov_no0_x15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x15(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x15(sa, ia, sc, ic, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x15(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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

integer :: s_o2, i_o2, s_o3, i_o3, s_o1, i_o1, s_i, i_i, s_k, i_k
! S2(i,k,a,c)  <-- 
! (    1.00000000) T2(o2,o3,o1,c) X(i,k,o3,o2,o1,a) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_o1,s_c) .and. &
IEOR(IEOR(s_i,s_k),s_o3) == IEOR(IEOR(s_o2,s_o1),s_a)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
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

end subroutine g_sigma_oovv_ooov_no1_x15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x16(sc, ic, V2, X, nir, nsym, psym)

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

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_oovv_ooov_no0_x16(sc, ic, h2_i, xaaaaa, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x16(s_c, i_c, V2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o4, i_o4, s_k, i_k, s_o3, i_o3, s_o5, i_o5
integer :: s_o2, i_o2, s_o1, i_o1
! X(i,o4,k,o3,o1,c)  <-- 
! (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(c,o1,o2,o5) 
do s_i = 0, nir-1
do s_o4 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o5 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_o1),s_c) .and. & 
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_o5),s_o2) .and. &
IEOR(s_c,s_o1) == IEOR(s_o2,s_o5)) then
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

end subroutine g_sigma_oovv_ooov_no0_x16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x16(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x16(sa, ia, sc, ic, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x16(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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

integer :: s_o3, i_o3, s_o4, i_o4, s_o1, i_o1, s_i, i_i, s_k, i_k
! S2(i,k,a,c)  <-- 
! (    1.00000000) T2(o3,o4,o1,a) X(i,o4,k,o3,o1,c) 
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o3,s_o4) == IEOR(s_o1,s_a) .and. &
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_o1),s_c)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o4, s_o3)%array(i_o1, i_o4, i_o3) & 
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

end subroutine g_sigma_oovv_ooov_no1_x16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x17(sa, ia, V2, X, nir, nsym, psym)

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

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_oovv_ooov_no0_x17(sa, ia, h2_i, xaaaaa, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x17(s_a, i_a, V2_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o4, i_o4, s_k, i_k, s_o3, i_o3, s_o5, i_o5
integer :: s_o2, i_o2, s_o1, i_o1
! X(i,o4,k,o3,o1,a)  <-- 
! (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(a,o1,o2,o5) 
do s_i = 0, nir-1
do s_o4 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o5 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_o1),s_a) .and. & 
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_o5),s_o2) .and. &
IEOR(s_a,s_o1) == IEOR(s_o2,s_o5)) then
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

end subroutine g_sigma_oovv_ooov_no0_x17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x17(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x17(sa, ia, sc, ic, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x17(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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

integer :: s_o4, i_o4, s_o3, i_o3, s_o1, i_o1, s_i, i_i, s_k, i_k
! S2(i,k,a,c)  <-- 
! (    1.00000000) T2(o4,o3,o1,c) X(i,o4,k,o3,o1,a) 
do s_o4 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o4,s_o3) == IEOR(s_o1,s_c) .and. &
IEOR(IEOR(s_i,s_o4),s_k) == IEOR(IEOR(s_o3,s_o1),s_a)) then
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
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

end subroutine g_sigma_oovv_ooov_no1_x17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x18(sc, ic, sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc

call set_symblock_xaaaav(sleft, x, nir, nsym, psym) ! -> xaaaav (allocate) 
call g_sigma_oovv_ooov_no0_x18(sc, ic, sv1, iv1, av2_i, h2_i, xaaaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaav)

end subroutine g_if_sigma_oovv_ooov_no0_x18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x18(s_c, i_c, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_o4, i_o4, s_o3, i_o3, s_a, i_a
! X(o1,o2,o4,o3,c,a)  <-- 
! (    1.00000000)  T2(o1,o2,o4,v1) V2(c,v1,o3,a) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o3 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(IEOR(s_o1,s_o2),s_o4) == IEOR(IEOR(s_o3,s_c),s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_o4,s_v1) .and. &
IEOR(s_c,s_v1) == IEOR(s_o3,s_a)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o3, s_o4, s_o2, s_o1)%array(i_a, i_o3, i_o4, i_o2, i_o1) =  &
    X_(s_a, s_o3, s_o4, s_o2, s_o1)%array(i_a, i_o3, i_o4, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_o4, s_o2, s_o1)%array(i_o4, i_o2, i_o1) & 
  * V2_(s_a, s_o3, s_v1)%array(i_a, i_o3, i_v1)
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

end subroutine g_sigma_oovv_ooov_no0_x18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x18(sc, ic, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sc

call set_symblock_xaaaav(sleft, x, nir, nsym, psym) ! -> xaaaav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x18(sc, ic, xaaaav, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x18(s_c, i_c, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1, s_a, i_a
! S2(i,k,a,c)  <-- 
! (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o1,o2,o4,o3,c,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(IEOR(s_o1,s_o2),s_o4) == IEOR(IEOR(s_o3,s_c),s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * X_(s_a, s_o3, s_o4, s_o2, s_o1)%array(i_a, i_o3, i_o4, i_o2, i_o1)
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

end subroutine g_sigma_oovv_ooov_no1_x18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x19(sc, ic, sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc

call set_symblock_xaaaav(sleft, x, nir, nsym, psym) ! -> xaaaav (allocate) 
call g_sigma_oovv_ooov_no0_x19(sc, ic, sv1, iv1, av2_i, h2_i, xaaaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaav)

end subroutine g_if_sigma_oovv_ooov_no0_x19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x19(s_c, i_c, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o4, i_o4, s_o2, i_o2, s_o1, i_o1, s_o3, i_o3, s_a, i_a
! X(o4,o2,o1,o3,c,a)  <-- 
! (    1.00000000)  T2(o4,o2,o1,v1) V2(c,o3,a,v1) 
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(IEOR(s_o4,s_o2),s_o1) == IEOR(IEOR(s_o3,s_c),s_a) .and. & 
IEOR(s_o4,s_o2) == IEOR(s_o1,s_v1) .and. &
IEOR(s_c,s_o3) == IEOR(s_a,s_v1)) then
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o3, s_o1, s_o2, s_o4)%array(i_a, i_o3, i_o1, i_o2, i_o4) =  &
    X_(s_a, s_o3, s_o1, s_o2, s_o4)%array(i_a, i_o3, i_o1, i_o2, i_o4) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o2, s_o4)%array(i_o1, i_o2, i_o4) & 
  * V2_(s_v1, s_a, s_o3)%array(i_v1, i_a, i_o3)
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

end subroutine g_sigma_oovv_ooov_no0_x19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x19(sc, ic, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sc

call set_symblock_xaaaav(sleft, x, nir, nsym, psym) ! -> xaaaav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x19(sc, ic, xaaaav, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x19(s_c, i_c, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o3, i_o3, s_o1, i_o1
integer :: s_o4, i_o4, s_a, i_a
! S2(i,k,a,c)  <-- 
! (    1.00000000) D3(i,o2,k,o3,o1,o4) X(o4,o2,o1,o3,c,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_i,s_o2),s_k) == IEOR(IEOR(s_o3,s_o1),s_o4) .and. &
IEOR(IEOR(s_o4,s_o2),s_o1) == IEOR(IEOR(s_o3,s_c),s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o4, s_o1, s_o3, s_k, s_o2, s_i)%array(i_o4, i_o1, i_o3, i_k, i_o2, i_i) & 
  * X_(s_a, s_o3, s_o1, s_o2, s_o4)%array(i_a, i_o3, i_o1, i_o2, i_o4)
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

end subroutine g_sigma_oovv_ooov_no1_x19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x20(sc, ic, sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_ooov_no0_x20(sc, ic, sv1, iv1, av2_i, h2_i, xaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaav)

end subroutine g_if_sigma_oovv_ooov_no0_x20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x20(s_c, i_c, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_o1, i_o1, s_a, i_a
! X(o3,o2,c,a)  <-- 
! (    1.00000000)  T2(o3,o2,o1,v1) V2(c,v1,o1,a) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o3,s_o2) == IEOR(s_c,s_a) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_o1,s_v1) .and. &
IEOR(s_c,s_v1) == IEOR(s_o1,s_a)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o2, s_o3)%array(i_a, i_o2, i_o3) =  &
    X_(s_a, s_o2, s_o3)%array(i_a, i_o2, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o2, s_o3)%array(i_o1, i_o2, i_o3) & 
  * V2_(s_a, s_o1, s_v1)%array(i_a, i_o1, i_v1)
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

end subroutine g_sigma_oovv_ooov_no0_x20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x20(sc, ic, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x20(sc, ic, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x20(s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_a, i_a
! S2(i,k,a,c)  <-- 
! (    1.00000000) D2(i,o3,k,o2) X(o3,o2,c,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. &
IEOR(s_o3,s_o2) == IEOR(s_c,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) & 
  * X_(s_a, s_o2, s_o3)%array(i_a, i_o2, i_o3)
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

end subroutine g_sigma_oovv_ooov_no1_x20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x21(sc, ic, sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_ooov_no0_x21(sc, ic, sv1, iv1, av2_i, h2_i, xaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaav)

end subroutine g_if_sigma_oovv_ooov_no0_x21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no0_x21(s_c, i_c, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_o3, i_o3, s_a, i_a
! X(o2,o1,c,a)  <-- 
! (    1.00000000)  T2(o2,o1,o3,v1) V2(c,o3,a,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_c,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_o3,s_v1) .and. &
IEOR(s_c,s_o3) == IEOR(s_a,s_v1)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) =  &
    X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_o3, s_o1, s_o2)%array(i_o3, i_o1, i_o2) & 
  * V2_(s_v1, s_a, s_o3)%array(i_v1, i_a, i_o3)
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

end subroutine g_sigma_oovv_ooov_no0_x21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x21(sc, ic, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x21(sc, ic, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ooov_no1_x21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ooov_no1_x21(s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2, s_a, i_a
! S2(i,k,a,c)  <-- 
! (    1.00000000) D2(i,o1,k,o2) X(o2,o1,c,a) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_o2) .and. &
IEOR(s_o2,s_o1) == IEOR(s_c,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i) & 
  * X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2)
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

end subroutine g_sigma_oovv_ooov_no1_x21

