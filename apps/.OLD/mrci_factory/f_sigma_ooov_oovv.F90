#include "../f_ct.fh"


!       #                #########      #     #   # 
!  ########## ##########         #   #######  #   # 
!      #    #         #          #    # #     #   # 
!      #    #        #   ########     # #     #   # 
!     #     #     # #           #  ##########    #  
!    #   # #       #            #       #       #   
!   #     #         #    ########       #     ##    



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x0(sa, ia, T2, h, X, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no0_x0(sa, ia, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x0(s_a, i_a, T2_, h_, X_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1, s_o3, i_o3
! X(o1,o2,o3,a)  <-- 
! (    1.00000000)  T2(o1,o2,v1,a) h(o3,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_v1,s_a) .and. &
IEOR(s_o3,s_v1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1) =  &
    X_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1) & 
  * h_(s_v1, s_o3)%array(i_v1, i_o3)
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

end subroutine g_sigma_ooov_oovv_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x0(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x0(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x0(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

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
! (    1.00000000) D3(i,m,k,o2,o3,o1) X(o1,o2,o3,a) 
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
  + 1.00000000d+00 & 
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

end subroutine g_sigma_ooov_oovv_no1_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x1(sa, ia, T2, h, X, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no0_x1(sa, ia, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x1(s_a, i_a, T2_, h_, X_, nir, nsym, psym)

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

integer :: s_o2, i_o2, s_o1, i_o1, s_v1, i_v1, s_m, i_m
! X(o2,o1,m,a)  <-- 
! (    1.00000000)  T2(o2,o1,v1,a) h(m,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_v1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_m,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_v1,s_a) .and. &
IEOR(s_m,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) =  &
    X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o1, s_o2)%array(i_v1, i_o1, i_o2) & 
  * h_(s_v1, s_m)%array(i_v1, i_m)
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

end subroutine g_sigma_ooov_oovv_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x1(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x1(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x1(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

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
! (    1.00000000) D2(i,o2,k,o1) X(o2,o1,m,a) 
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
  + 1.00000000d+00 & 
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

end subroutine g_sigma_ooov_oovv_no1_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x2(sa, ia, sv1, iv1, T2, h, X, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no0_x2(sa, ia, sv1, iv1, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x2(s_a, i_a, s_v1, i_v1, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o3, i_o3, s_o1, i_o1
! X(o2,o3,o1,a)  <-- 
! (    1.00000000)  T2(o2,o3,a,v1) h(o1,v1) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_o2,s_o3) == IEOR(s_o1,s_a) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_a,s_v1) .and. &
IEOR(s_o1,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2) =  &
    X_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o3, s_o2)%array(i_a, i_o3, i_o2) & 
  * h_(s_v1, s_o1)%array(i_v1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_oovv_no0_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x2(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x2(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x2(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

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

end subroutine g_sigma_ooov_oovv_no1_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x3(sa, ia, sv1, iv1, T2, h, X, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no0_x3(sa, ia, sv1, iv1, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x3(s_a, i_a, s_v1, i_v1, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_m, i_m
! X(o2,o1,m,a)  <-- 
! (    1.00000000)  T2(o2,o1,a,v1) h(m,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_m,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_a,s_v1) .and. &
IEOR(s_m,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) =  &
    X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) & 
  * h_(s_v1, s_m)%array(i_v1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_oovv_no0_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x3(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x3(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x3(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

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

end subroutine g_sigma_ooov_oovv_no1_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_y4(sc1, ic1, V2, Y0, nir, nsym, psym)

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
call g_sigma_ooov_oovv_y4(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_ooov_oovv_y4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_y4(s_c1, i_c1, V2_, Y0_, nir, nsym, psym)

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

integer :: s_o3, i_o3, s_v1, i_v1
! Y0(o3,v1)  <-- 
! (    1.00000000)  V2(c1,c1,o3,v1) 
do s_o3 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_o3,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o3,s_v1)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y0_(s_v1, s_o3)%array(i_v1, i_o3) =  &
    Y0_(s_v1, s_o3)%array(i_v1, i_o3) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_o3, s_c1)%array(i_v1, i_o3, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_oovv_y4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x4(sa, ia, T2, Y0, X, nir, nsym, psym)

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
call set_symblock_Yav(sleft2, Y0, nir, nsym, psym) ! -> Yav (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_oovv_no0_x4(sa, ia, av2_i, Yav, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yav)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x4(s_a, i_a, T2_, Y0_, X_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1, s_o3, i_o3
! X(o1,o2,o3,a)  <-- 
! (    1.00000000)  T2(o1,o2,v1,a) Y0(o3,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_v1,s_a) .and. &
IEOR(s_o3,s_v1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1) =  &
    X_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1) & 
  * Y0_(s_v1, s_o3)%array(i_v1, i_o3)
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

end subroutine g_sigma_ooov_oovv_no0_x4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x4(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x4(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x4(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

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

end subroutine g_sigma_ooov_oovv_no1_x4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_y5(sc1, ic1, V2, Y1, nir, nsym, psym)

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
call g_sigma_ooov_oovv_y5(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_ooov_oovv_y5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_y5(s_c1, i_c1, V2_, Y1_, nir, nsym, psym)

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

integer :: s_o3, i_o3, s_v1, i_v1
! Y1(o3,v1)  <-- 
! (    1.00000000)  V2(c1,o3,c1,v1) 
do s_o3 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_o3,s_v1) == 0 .and. & 
IEOR(s_c1,s_o3) == IEOR(s_c1,s_v1)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y1_(s_v1, s_o3)%array(i_v1, i_o3) =  &
    Y1_(s_v1, s_o3)%array(i_v1, i_o3) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_o3)%array(i_v1, i_c1, i_o3)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_oovv_y5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x5(sa, ia, T2, Y1, X, nir, nsym, psym)

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
call set_symblock_Yav(sleft2, Y1, nir, nsym, psym) ! -> Yav (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_oovv_no0_x5(sa, ia, av2_i, Yav, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yav)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x5(s_a, i_a, T2_, Y1_, X_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1, s_o3, i_o3
! X(o1,o2,o3,a)  <-- 
! (    1.00000000)  T2(o1,o2,v1,a) Y1(o3,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_v1,s_a) .and. &
IEOR(s_o3,s_v1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1) =  &
    X_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1) & 
  * Y1_(s_v1, s_o3)%array(i_v1, i_o3)
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

end subroutine g_sigma_ooov_oovv_no0_x5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x5(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x5(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x5(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

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

end subroutine g_sigma_ooov_oovv_no1_x5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_y6(sc1, ic1, V2, Y2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_y6(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_ooov_oovv_y6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_y6(s_c1, i_c1, V2_, Y2_, nir, nsym, psym)

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

integer :: s_m, i_m, s_v1, i_v1
! Y2(m,v1)  <-- 
! (    1.00000000)  V2(c1,c1,m,v1) 
do s_m = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_m,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_m,s_v1)) then
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y2_(s_v1, s_m)%array(i_v1, i_m) =  &
    Y2_(s_v1, s_m)%array(i_v1, i_m) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_m, s_c1)%array(i_v1, i_m, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_oovv_y6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x6(sa, ia, T2, Y2, X, nir, nsym, psym)

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
call set_symblock_Yav(sleft2, Y2, nir, nsym, psym) ! -> Yav (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_oovv_no0_x6(sa, ia, av2_i, Yav, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yav)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x6(s_a, i_a, T2_, Y2_, X_, nir, nsym, psym)

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

integer :: s_o2, i_o2, s_o1, i_o1, s_v1, i_v1, s_m, i_m
! X(o2,o1,m,a)  <-- 
! (    1.00000000)  T2(o2,o1,v1,a) Y2(m,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_v1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_m,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_v1,s_a) .and. &
IEOR(s_m,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) =  &
    X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o1, s_o2)%array(i_v1, i_o1, i_o2) & 
  * Y2_(s_v1, s_m)%array(i_v1, i_m)
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

end subroutine g_sigma_ooov_oovv_no0_x6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x6(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x6(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x6(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

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

end subroutine g_sigma_ooov_oovv_no1_x6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_y7(sc1, ic1, V2, Y3, nir, nsym, psym)

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
call g_sigma_ooov_oovv_y7(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_ooov_oovv_y7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_y7(s_c1, i_c1, V2_, Y3_, nir, nsym, psym)

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

integer :: s_m, i_m, s_v1, i_v1
! Y3(m,v1)  <-- 
! (    1.00000000)  V2(c1,m,c1,v1) 
do s_m = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_m,s_v1) == 0 .and. & 
IEOR(s_c1,s_m) == IEOR(s_c1,s_v1)) then
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y3_(s_v1, s_m)%array(i_v1, i_m) =  &
    Y3_(s_v1, s_m)%array(i_v1, i_m) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_m)%array(i_v1, i_c1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_oovv_y7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x7(sa, ia, T2, Y3, X, nir, nsym, psym)

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
call set_symblock_Yav(sleft2, Y3, nir, nsym, psym) ! -> Yav (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_oovv_no0_x7(sa, ia, av2_i, Yav, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yav)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x7(s_a, i_a, T2_, Y3_, X_, nir, nsym, psym)

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

integer :: s_o2, i_o2, s_o1, i_o1, s_v1, i_v1, s_m, i_m
! X(o2,o1,m,a)  <-- 
! (    1.00000000)  T2(o2,o1,v1,a) Y3(m,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_v1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_m,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_v1,s_a) .and. &
IEOR(s_m,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) =  &
    X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o1, s_o2)%array(i_v1, i_o1, i_o2) & 
  * Y3_(s_v1, s_m)%array(i_v1, i_m)
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

end subroutine g_sigma_ooov_oovv_no0_x7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x7(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x7(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x7(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

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

end subroutine g_sigma_ooov_oovv_no1_x7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_y8(sc1, ic1, V2, Y4, nir, nsym, psym)

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
call g_sigma_ooov_oovv_y8(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_ooov_oovv_y8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_y8(s_c1, i_c1, V2_, Y4_, nir, nsym, psym)

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

integer :: s_o3, i_o3, s_v1, i_v1
! Y4(o3,v1)  <-- 
! (    1.00000000)  V2(c1,c1,o3,v1) 
do s_o3 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_o3,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o3,s_v1)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y4_(s_v1, s_o3)%array(i_v1, i_o3) =  &
    Y4_(s_v1, s_o3)%array(i_v1, i_o3) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_o3, s_c1)%array(i_v1, i_o3, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_oovv_y8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x8(sa, ia, sv1, iv1, T2, Y4, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y4(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yav(sleft2, Y4, nir, nsym, psym) ! -> Yav (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_oovv_no0_x8(sa, ia, sv1, iv1, av2_i, Yav, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yav)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x8(s_a, i_a, s_v1, i_v1, T2_, Y4_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y4_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_o3, i_o3
! X(o2,o1,o3,a)  <-- 
! (    1.00000000)  T2(o2,o1,a,v1) Y4(o3,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_a,s_v1) .and. &
IEOR(s_o3,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o1, s_o2)%array(i_o3, i_o1, i_o2) =  &
    X_(s_o3, s_o1, s_o2)%array(i_o3, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) & 
  * Y4_(s_v1, s_o3)%array(i_v1, i_o3)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_oovv_no0_x8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x8(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x8(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x8(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

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
! (    2.00000000) D3(i,m,k,o2,o3,o1) X(o2,o1,o3,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o2,s_o3),s_o1) .and. &
IEOR(s_o2,s_o1) == IEOR(s_o3,s_a)) then
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
  * X_(s_o3, s_o1, s_o2)%array(i_o3, i_o1, i_o2)
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

end subroutine g_sigma_ooov_oovv_no1_x8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_y9(sc1, ic1, V2, Y5, nir, nsym, psym)

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
call g_sigma_ooov_oovv_y9(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_ooov_oovv_y9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_y9(s_c1, i_c1, V2_, Y5_, nir, nsym, psym)

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

integer :: s_o3, i_o3, s_v1, i_v1
! Y5(o3,v1)  <-- 
! (    1.00000000)  V2(c1,o3,c1,v1) 
do s_o3 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_o3,s_v1) == 0 .and. & 
IEOR(s_c1,s_o3) == IEOR(s_c1,s_v1)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y5_(s_v1, s_o3)%array(i_v1, i_o3) =  &
    Y5_(s_v1, s_o3)%array(i_v1, i_o3) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_o3)%array(i_v1, i_c1, i_o3)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_oovv_y9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x9(sa, ia, sv1, iv1, T2, Y5, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y5(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yav(sleft2, Y5, nir, nsym, psym) ! -> Yav (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_oovv_no0_x9(sa, ia, sv1, iv1, av2_i, Yav, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yav)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x9(s_a, i_a, s_v1, i_v1, T2_, Y5_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y5_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_o3, i_o3
! X(o2,o1,o3,a)  <-- 
! (    1.00000000)  T2(o2,o1,a,v1) Y5(o3,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_a,s_v1) .and. &
IEOR(s_o3,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o1, s_o2)%array(i_o3, i_o1, i_o2) =  &
    X_(s_o3, s_o1, s_o2)%array(i_o3, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) & 
  * Y5_(s_v1, s_o3)%array(i_v1, i_o3)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_oovv_no0_x9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x9(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x9(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x9(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

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
! (   -1.00000000) D3(i,m,k,o2,o3,o1) X(o2,o1,o3,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o2,s_o3),s_o1) .and. &
IEOR(s_o2,s_o1) == IEOR(s_o3,s_a)) then
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
  * X_(s_o3, s_o1, s_o2)%array(i_o3, i_o1, i_o2)
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

end subroutine g_sigma_ooov_oovv_no1_x9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_y10(sc1, ic1, V2, Y6, nir, nsym, psym)

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
call g_sigma_ooov_oovv_y10(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_ooov_oovv_y10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_y10(s_c1, i_c1, V2_, Y6_, nir, nsym, psym)

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

integer :: s_m, i_m, s_v1, i_v1
! Y6(m,v1)  <-- 
! (    1.00000000)  V2(c1,c1,m,v1) 
do s_m = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_m,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_m,s_v1)) then
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y6_(s_v1, s_m)%array(i_v1, i_m) =  &
    Y6_(s_v1, s_m)%array(i_v1, i_m) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_m, s_c1)%array(i_v1, i_m, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_oovv_y10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x10(sa, ia, sv1, iv1, T2, Y6, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y6(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yav(sleft2, Y6, nir, nsym, psym) ! -> Yav (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_oovv_no0_x10(sa, ia, sv1, iv1, av2_i, Yav, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yav)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x10(s_a, i_a, s_v1, i_v1, T2_, Y6_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y6_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_m, i_m
! X(o1,o2,m,a)  <-- 
! (    1.00000000)  T2(o1,o2,a,v1) Y6(m,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_m,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_a,s_v1) .and. &
IEOR(s_m,s_v1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o2, s_o1)%array(i_m, i_o2, i_o1) =  &
    X_(s_m, s_o2, s_o1)%array(i_m, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1) & 
  * Y6_(s_v1, s_m)%array(i_v1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_oovv_no0_x10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x10(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x10(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x10(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

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
! (    2.00000000) D2(i,o2,k,o1) X(o1,o2,m,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_m,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 2.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
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

end subroutine g_sigma_ooov_oovv_no1_x10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_y11(sc1, ic1, V2, Y7, nir, nsym, psym)

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
call g_sigma_ooov_oovv_y11(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_ooov_oovv_y11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_y11(s_c1, i_c1, V2_, Y7_, nir, nsym, psym)

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

integer :: s_m, i_m, s_v1, i_v1
! Y7(m,v1)  <-- 
! (    1.00000000)  V2(c1,m,c1,v1) 
do s_m = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_m,s_v1) == 0 .and. & 
IEOR(s_c1,s_m) == IEOR(s_c1,s_v1)) then
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y7_(s_v1, s_m)%array(i_v1, i_m) =  &
    Y7_(s_v1, s_m)%array(i_v1, i_m) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_m)%array(i_v1, i_c1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_oovv_y11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x11(sa, ia, sv1, iv1, T2, Y7, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y7(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yav(sleft2, Y7, nir, nsym, psym) ! -> Yav (allocate) 
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_oovv_no0_x11(sa, ia, sv1, iv1, av2_i, Yav, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yav)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x11(s_a, i_a, s_v1, i_v1, T2_, Y7_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y7_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_m, i_m
! X(o1,o2,m,a)  <-- 
! (    1.00000000)  T2(o1,o2,a,v1) Y7(m,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_m,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_a,s_v1) .and. &
IEOR(s_m,s_v1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o2, s_o1)%array(i_m, i_o2, i_o1) =  &
    X_(s_m, s_o2, s_o1)%array(i_m, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1) & 
  * Y7_(s_v1, s_m)%array(i_v1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_oovv_no0_x11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x11(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x11(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x11(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

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
! (   -1.00000000) D2(i,o2,k,o1) X(o1,o2,m,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_m,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
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

end subroutine g_sigma_ooov_oovv_no1_x11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x12(sa, ia, so2, io2, so5, io5, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, so2, io2, so5, io5
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(so2, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(so2,IEOR(so5,sa))

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_oovv_no0_x12(sa, ia, so2, io2, so5, io5, av2_i, h2_i, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x12(s_a, i_a, s_o2, i_o2, s_o5, i_o5, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_o2, s_o2
integer, intent(in) :: i_o5, s_o5
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3, s_v1, i_v1, s_o4, i_o4
! X(o1,o3,o2,o5,o4,a)  <-- 
! (    1.00000000)  T2(o1,o3,v1,a) V2(o2,o5,o4,v1) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_v1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(IEOR(s_o1,s_o3),s_o2) == IEOR(IEOR(s_o5,s_o4),s_a) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_v1,s_a) .and. &
IEOR(s_o2,s_o5) == IEOR(s_o4,s_v1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o3, s_o1)%array(i_o4, i_o3, i_o1) =  &
    X_(s_o4, s_o3, s_o1)%array(i_o4, i_o3, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o3, s_o1)%array(i_v1, i_o3, i_o1) & 
  * V2_(s_v1, s_o4, s_o5)%array(i_v1, i_o4, i_o5)
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

end subroutine g_sigma_ooov_oovv_no0_x12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x12(sa, ia, so2, io2, so5, io5, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, so2, io2, so5, io5
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = IEOR(so2,IEOR(so5,sa))

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_oovv_no1_x12(sa, ia, so2, io2, so5, io5, xaaa, av2_i2, d4_ij, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x12(s_a, i_a, s_o2, i_o2, s_o5, i_o5, X_, S2_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o2, s_o2
integer, intent(in) :: i_o5, s_o5
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_m, i_m, s_i, i_i, s_o3, i_o3, s_k, i_k, s_o1, i_o1
integer :: s_o4, i_o4
! S2(i,k,m,a)  <-- 
! (    1.00000000) D4(o2,o5,m,i,o3,k,o1,o4) X(o1,o3,o2,o5,o4,a) 
do s_m = 0, nir-1
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_o2,s_o5),IEOR(s_m,s_i)) == IEOR(IEOR(s_o3,s_k),IEOR(s_o1,s_o4)) .and. &
IEOR(IEOR(s_o1,s_o3),s_o2) == IEOR(IEOR(s_o5,s_o4),s_a)) then
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D4_(s_o4, s_o1, s_k, s_o3, s_i, s_m)%array(i_o4, i_o1, i_k, i_o3, i_i, i_m) & 
  * X_(s_o4, s_o3, s_o1)%array(i_o4, i_o3, i_o1)
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

end subroutine g_sigma_ooov_oovv_no1_x12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x13(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

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

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_oovv_no0_x13(sa, ia, sv1, iv1, av2_i, h2_i, xaaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x13(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

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
! (    1.00000000)  T2(o1,o2,v1,a) V2(v1,o4,m,o3) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_m = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(IEOR(s_o1,s_o2),s_o4) == IEOR(IEOR(s_m,s_o3),s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_v1,s_a) .and. &
IEOR(s_v1,s_o4) == IEOR(s_m,s_o3)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_m, s_o4, s_o2, s_o1)%array(i_o3, i_m, i_o4, i_o2, i_o1) =  &
    X_(s_o3, s_m, s_o4, s_o2, s_o1)%array(i_o3, i_m, i_o4, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1) & 
  * V2_(s_o3, s_m, s_o4)%array(i_o3, i_m, i_o4)
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

end subroutine g_sigma_ooov_oovv_no0_x13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x13(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x13(sa, ia, xaaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x13(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

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

end subroutine g_sigma_ooov_oovv_no1_x13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x14(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

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

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_oovv_no0_x14(sa, ia, sv1, iv1, av2_i, h2_i, xaaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x14(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

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
! (    1.00000000)  T2(o3,o2,v1,a) V2(v1,m,o1,o4) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(IEOR(s_o3,s_o2),s_m) == IEOR(IEOR(s_o1,s_o4),s_a) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_v1,s_a) .and. &
IEOR(s_v1,s_m) == IEOR(s_o1,s_o4)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o1, s_m, s_o2, s_o3)%array(i_o4, i_o1, i_m, i_o2, i_o3) =  &
    X_(s_o4, s_o1, s_m, s_o2, s_o3)%array(i_o4, i_o1, i_m, i_o2, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o2, s_o3)%array(i_v1, i_o2, i_o3) & 
  * V2_(s_o4, s_o1, s_m)%array(i_o4, i_o1, i_m)
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

end subroutine g_sigma_ooov_oovv_no0_x14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x14(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x14(sa, ia, xaaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x14(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

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

end subroutine g_sigma_ooov_oovv_no1_x14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x15(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

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

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_oovv_no0_x15(sa, ia, sv1, iv1, av2_i, h2_i, xaaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x15(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o5, i_o5, s_o2, i_o2, s_o1, i_o1, s_o4, i_o4
! X(o3,o5,o2,o1,o4,a)  <-- 
! (    1.00000000)  T2(o3,o5,a,v1) V2(v1,o2,o1,o4) 
do s_o3 = 0, nir-1
do s_o5 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(IEOR(s_o3,s_o5),s_o2) == IEOR(IEOR(s_o1,s_o4),s_a) .and. & 
IEOR(s_o3,s_o5) == IEOR(s_a,s_v1) .and. &
IEOR(s_v1,s_o2) == IEOR(s_o1,s_o4)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o1, s_o2, s_o5, s_o3)%array(i_o4, i_o1, i_o2, i_o5, i_o3) =  &
    X_(s_o4, s_o1, s_o2, s_o5, s_o3)%array(i_o4, i_o1, i_o2, i_o5, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o5, s_o3)%array(i_a, i_o5, i_o3) & 
  * V2_(s_o4, s_o1, s_o2)%array(i_o4, i_o1, i_o2)
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

end subroutine g_sigma_ooov_oovv_no0_x15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x15(sa, ia, si, ii, sm, im, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x15(sa, ia, si, ii, sm, im, xaaaaa, av2_i2, d4_ij, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x15(s_a, i_a, s_i, i_i, s_m, i_m, X_, S2_, D4_, nir, nsym, psym)

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
! (    1.00000000) D4(i,m,k,o3,o1,o4,o2,o5) X(o3,o5,o2,o1,o4,a) 
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
do s_o5 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),IEOR(s_k,s_o3)) == IEOR(IEOR(s_o1,s_o4),IEOR(s_o2,s_o5)) .and. &
IEOR(IEOR(s_o3,s_o5),s_o2) == IEOR(IEOR(s_o1,s_o4),s_a)) then
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
  * X_(s_o4, s_o1, s_o2, s_o5, s_o3)%array(i_o4, i_o1, i_o2, i_o5, i_o3)
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

end subroutine g_sigma_ooov_oovv_no1_x15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x16(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

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

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_oovv_no0_x16(sa, ia, sv1, iv1, av2_i, h2_i, xaaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x16(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_o4, i_o4, s_m, i_m, s_o3, i_o3
! X(o2,o1,o4,m,o3,a)  <-- 
! (    1.00000000)  T2(o2,o1,a,v1) V2(v1,o4,m,o3) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_m = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(IEOR(s_o2,s_o1),s_o4) == IEOR(IEOR(s_m,s_o3),s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_a,s_v1) .and. &
IEOR(s_v1,s_o4) == IEOR(s_m,s_o3)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_m, s_o4, s_o1, s_o2)%array(i_o3, i_m, i_o4, i_o1, i_o2) =  &
    X_(s_o3, s_m, s_o4, s_o1, s_o2)%array(i_o3, i_m, i_o4, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) & 
  * V2_(s_o3, s_m, s_o4)%array(i_o3, i_m, i_o4)
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

end subroutine g_sigma_ooov_oovv_no0_x16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x16(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x16(sa, ia, xaaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x16(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

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
! (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o2,o1,o4,m,o3,a) 
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
IEOR(IEOR(s_o2,s_o1),s_o4) == IEOR(IEOR(s_m,s_o3),s_a)) then
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
  * X_(s_o3, s_m, s_o4, s_o1, s_o2)%array(i_o3, i_m, i_o4, i_o1, i_o2)
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

end subroutine g_sigma_ooov_oovv_no1_x16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x17(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

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

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_oovv_no0_x17(sa, ia, sv1, iv1, av2_i, h2_i, xaaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x17(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o3, i_o3, s_m, i_m, s_o1, i_o1, s_o4, i_o4
! X(o2,o3,m,o1,o4,a)  <-- 
! (    1.00000000)  T2(o2,o3,a,v1) V2(v1,m,o1,o4) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_m = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(IEOR(s_o2,s_o3),s_m) == IEOR(IEOR(s_o1,s_o4),s_a) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_a,s_v1) .and. &
IEOR(s_v1,s_m) == IEOR(s_o1,s_o4)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o1, s_m, s_o3, s_o2)%array(i_o4, i_o1, i_m, i_o3, i_o2) =  &
    X_(s_o4, s_o1, s_m, s_o3, s_o2)%array(i_o4, i_o1, i_m, i_o3, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o3, s_o2)%array(i_a, i_o3, i_o2) & 
  * V2_(s_o4, s_o1, s_m)%array(i_o4, i_o1, i_m)
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

end subroutine g_sigma_ooov_oovv_no0_x17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x17(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x17(sa, ia, xaaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x17(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

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
! (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o2,o3,m,o1,o4,a) 
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
IEOR(IEOR(s_o2,s_o3),s_m) == IEOR(IEOR(s_o1,s_o4),s_a)) then
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
  * X_(s_o4, s_o1, s_m, s_o3, s_o2)%array(i_o4, i_o1, i_m, i_o3, i_o2)
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

end subroutine g_sigma_ooov_oovv_no1_x17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x18(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no0_x18(sa, ia, sv1, iv1, av2_i, h2_i, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x18(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o2, i_o2, s_v2, i_v2, s_o3, i_o3
! X(o1,o2,o3,a)  <-- 
! (    1.00000000)  T2(o1,o2,v2,v1) V2(a,v1,o3,v2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v2 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_o3,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_v2,s_v1) .and. &
IEOR(s_a,s_v1) == IEOR(s_o3,s_v2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v2 = psym(I_BEGIN, I_V, s_v2), psym(I_END, I_V, s_v2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1) =  &
    X_(s_o3, s_o2, s_o1)%array(i_o3, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_v2, s_o2, s_o1)%array(i_v2, i_o2, i_o1) & 
  * V2_(s_v2, s_o3, s_v1)%array(i_v2, i_o3, i_v1)
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

end subroutine g_sigma_ooov_oovv_no0_x18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x18(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x18(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x18(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

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
! (    1.00000000) D3(i,m,k,o2,o3,o1) X(o1,o2,o3,a) 
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
  + 1.00000000d+00 & 
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

end subroutine g_sigma_ooov_oovv_no1_x18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x19(sa, ia, sv1, iv1, T2, V2, X, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no0_x19(sa, ia, sv1, iv1, av2_i, h2_i, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x19(s_a, i_a, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

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

integer :: s_o2, i_o2, s_o1, i_o1, s_v2, i_v2, s_m, i_m
! X(o2,o1,m,a)  <-- 
! (    1.00000000)  T2(o2,o1,v2,v1) V2(a,v1,m,v2) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_v2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_m,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_v2,s_v1) .and. &
IEOR(s_a,s_v1) == IEOR(s_m,s_v2)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_v2 = psym(I_BEGIN, I_V, s_v2), psym(I_END, I_V, s_v2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) =  &
    X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_v2, s_o1, s_o2)%array(i_v2, i_o1, i_o2) & 
  * V2_(s_v2, s_m, s_v1)%array(i_v2, i_m, i_v1)
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

end subroutine g_sigma_ooov_oovv_no0_x19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x19(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x19(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x19(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

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
! (    1.00000000) D2(i,o2,k,o1) X(o2,o1,m,a) 
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
  + 1.00000000d+00 & 
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

end subroutine g_sigma_ooov_oovv_no1_x19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x20(sa, ia, sv2, iv2, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv2, iv2
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv2, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_oovv_no0_x20(sa, ia, sv2, iv2, av2_i, h2_i, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x20(s_a, i_a, s_v2, i_v2, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v2, s_v2
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o3, i_o3, s_v1, i_v1, s_o1, i_o1
! X(o2,o3,o1,a)  <-- 
! (    1.00000000)  T2(o2,o3,v1,v2) V2(a,v1,o1,v2) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_v1 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_o2,s_o3) == IEOR(s_o1,s_a) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_v1,s_v2) .and. &
IEOR(s_a,s_v1) == IEOR(s_o1,s_v2)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2) =  &
    X_(s_o1, s_o3, s_o2)%array(i_o1, i_o3, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o3, s_o2)%array(i_v1, i_o3, i_o2) & 
  * V2_(s_v2, s_o1, s_v1)%array(i_v2, i_o1, i_v1)
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

end subroutine g_sigma_ooov_oovv_no0_x20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x20(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x20(sa, ia, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x20(s_a, i_a, X_, S2_, D3_, nir, nsym, psym)

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

end subroutine g_sigma_ooov_oovv_no1_x20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no0_x21(sa, ia, sv2, iv2, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sv2, iv2
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv2, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_oovv_no0_x21(sa, ia, sv2, iv2, av2_i, h2_i, xaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_oovv_no0_x21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no0_x21(s_a, i_a, s_v2, i_v2, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v2, s_v2
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_v1, i_v1, s_m, i_m
! X(o2,o1,m,a)  <-- 
! (    1.00000000)  T2(o2,o1,v1,v2) V2(a,v1,m,v2) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_v1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_m,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_v1,s_v2) .and. &
IEOR(s_a,s_v1) == IEOR(s_m,s_v2)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) =  &
    X_(s_m, s_o1, s_o2)%array(i_m, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o1, s_o2)%array(i_v1, i_o1, i_o2) & 
  * V2_(s_v2, s_m, s_v1)%array(i_v2, i_m, i_v1)
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

end subroutine g_sigma_ooov_oovv_no0_x21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_oovv_no1_x21(sa, ia, X, S2, nir, nsym, psym)

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
call g_sigma_ooov_oovv_no1_x21(sa, ia, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_oovv_no1_x21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_oovv_no1_x21(s_a, i_a, X_, S2_, D2_, nir, nsym, psym)

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

end subroutine g_sigma_ooov_oovv_no1_x21

