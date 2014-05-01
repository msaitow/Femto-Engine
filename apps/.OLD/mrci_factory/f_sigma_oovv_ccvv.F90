#include "../f_ct.fh"


!        :::::::::: ::::::::::   :::   ::: ::::::::::: :::::::: 
!       :+:        :+:         :+:+: :+:+:    :+:    :+:    :+: 
!      +:+        +:+        +:+ +:+:+ +:+   +:+    +:+    +:+  
!     :#::+::#   +#++:++#   +#+  +:+  +#+   +#+    +#+    +:+   
!    +#+        +#+        +#+       +#+   +#+    +#+    +#+    
!   #+#        #+#        #+#       #+#   #+#    #+#    #+#     
!  ###        ########## ###       ###   ###     ########       



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ccvv_no0_x0(sc, ic, sc1, ic1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_ccvv_no0_x0(sc, ic, sc1, ic1, av2_i, h2_i, xaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaav)

end subroutine g_if_sigma_oovv_ccvv_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ccvv_no0_x0(s_c, i_c, s_c1, i_c1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_a, i_a, s_o1, i_o1, s_o2, i_o2
! X(o1,o2,a,c)  <-- 
! (    1.00000000)  T2(c2,c1,a,c) V2(c1,o1,c2,o2) 
do s_c2 = 0, nir-1
do s_a = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_a,s_c) .and. & 
IEOR(s_c2,s_c1) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_o1) == IEOR(s_c2,s_o2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1) =  &
    X_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_c1, s_c2)%array(i_a, i_c1, i_c2) & 
  * V2_(s_o2, s_c2, s_o1)%array(i_o2, i_c2, i_o1)
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

end subroutine g_sigma_oovv_ccvv_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ccvv_no1_x0(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_ccvv_no1_x0(sc, ic, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ccvv_no1_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ccvv_no1_x0(s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

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

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o1, i_o1, s_a, i_a
! S2(i,k,a,c)  <-- 
! (    1.00000000) D2(i,o2,k,o1) X(o1,o2,a,c) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_a,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
  * X_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
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

end subroutine g_sigma_oovv_ccvv_no1_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ccvv_no0_x1(sa, ia, sc, ic, sc1, ic1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sc,sa)

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call g_sigma_oovv_ccvv_no0_x1(sa, ia, sc, ic, sc1, ic1, av2_i, h2_i, xaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaa)

end subroutine g_if_sigma_oovv_ccvv_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ccvv_no0_x1(s_a, i_a, s_c, i_c, s_c1, i_c1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_o2, i_o2, s_o1, i_o1
! X(o2,o1,c,a)  <-- 
! (    1.00000000)  T2(c1,c2,c,a) V2(c1,o2,c2,o1) 
do s_c2 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_c,s_a) .and. & 
IEOR(s_c1,s_c2) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_o2) == IEOR(s_c2,s_o1)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_o2)%array(i_o1, i_o2) =  &
    X_(s_o1, s_o2)%array(i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_c, s_c2, s_c1)%array(i_c, i_c2, i_c1) & 
  * V2_(s_o1, s_c2, s_o2)%array(i_o1, i_c2, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ccvv_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ccvv_no1_x1(sa, ia, sc, ic, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = IEOR(sc,sa)

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ccvv_no1_x1(sa, ia, sc, ic, xaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_ccvv_no1_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_ccvv_no1_x1(s_a, i_a, s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2
! S2(i,k,a,c)  <-- 
! (    1.00000000) D2(i,o1,k,o2) X(o2,o1,c,a) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_o2) .and. &
IEOR(s_o2,s_o1) == IEOR(s_c,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i) & 
  * X_(s_o1, s_o2)%array(i_o1, i_o2)
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

end subroutine g_sigma_oovv_ccvv_no1_x1

