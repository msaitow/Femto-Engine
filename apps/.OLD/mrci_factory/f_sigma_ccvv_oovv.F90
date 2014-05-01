#include "../f_ct.fh"


!  8888888888                     888                  
!  888                            888                  
!  888                            888                  
!  8888888  .d88b.  88888b.d88b.  888888  .d88b.       
!  888     d8P  Y8b 888 "888 "88b 888    d88""88b  
!  888     88888888 888  888  888 888    888  888      
!  888     Y8b.     888  888  888 Y88b.  Y88..88P      
!  888      "Y8888  888  888  888  "Y888  "Y88P"   



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_oovv_no0_x0(sw, iw, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_xcaa(sleft, x, nir, nsym, psym) ! -> xcaa (allocate) 
call g_sigma_ccvv_oovv_no0_x0(sw, iw, h2_i, xcaa, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(xcaa)

end subroutine g_if_sigma_ccvv_oovv_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_oovv_no0_x0(s_w, i_w, V2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3, s_o2, i_o2, s_o4, i_o4, s_y, i_y
! X(w,y,o1,o2)  <-- 
! (    1.00000000)  D2(o1,o3,o2,o4) V2(w,o3,y,o4) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_o1,s_o2) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_w,s_o3) == IEOR(s_y,s_o4)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
X_(s_o2, s_o1, s_y)%array(i_o2, i_o1, i_y) =  &
    X_(s_o2, s_o1, s_y)%array(i_o2, i_o1, i_y) &
  + 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * V2_(s_o4, s_y, s_o3)%array(i_o4, i_y, i_o3)
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

end subroutine g_sigma_ccvv_oovv_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_oovv_no1_x0(sc, ic, sw, iw, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sw, iw
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_xcaa(sleft, x, nir, nsym, psym) ! -> xcaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_oovv_no1_x0(sc, ic, sw, iw, av2_i, xcaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xcaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_oovv_no1_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_oovv_no1_x0(s_c, i_c, s_w, i_w, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    1.00000000) T2(o1,o2,a,c) X(w,y,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_a,s_c) .and. &
IEOR(s_w,s_y) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1) & 
  * X_(s_o2, s_o1, s_y)%array(i_o2, i_o1, i_y)
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

end subroutine g_sigma_ccvv_oovv_no1_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_oovv_no0_x1(sw, iw, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_xcaa(sleft, x, nir, nsym, psym) ! -> xcaa (allocate) 
call g_sigma_ccvv_oovv_no0_x1(sw, iw, h2_i, xcaa, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(xcaa)

end subroutine g_if_sigma_ccvv_oovv_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_oovv_no0_x1(s_w, i_w, V2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3, s_o2, i_o2, s_o4, i_o4, s_y, i_y
! X(w,y,o3,o4)  <-- 
! (    1.00000000)  D2(o1,o3,o2,o4) V2(w,o2,y,o1) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_o3,s_o4) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_w,s_o2) == IEOR(s_y,s_o1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
X_(s_o4, s_o3, s_y)%array(i_o4, i_o3, i_y) =  &
    X_(s_o4, s_o3, s_y)%array(i_o4, i_o3, i_y) &
  + 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * V2_(s_o1, s_y, s_o2)%array(i_o1, i_y, i_o2)
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

end subroutine g_sigma_ccvv_oovv_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_oovv_no1_x1(sa, ia, sc, ic, sw, iw, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic, sw, iw
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_xcaa(sleft, x, nir, nsym, psym) ! -> xcaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_oovv_no1_x1(sa, ia, sc, ic, sw, iw, av2_i, xcaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xcaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_oovv_no1_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_oovv_no1_x1(s_a, i_a, s_c, i_c, s_w, i_w, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o4, i_o4, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    1.00000000) T2(o3,o4,c,a) X(w,y,o3,o4) 
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o3,s_o4) == IEOR(s_c,s_a) .and. &
IEOR(s_w,s_y) == IEOR(s_o3,s_o4)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * T2_(s_c, s_o4, s_o3)%array(i_c, i_o4, i_o3) & 
  * X_(s_o4, s_o3, s_y)%array(i_o4, i_o3, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_oovv_no1_x1

