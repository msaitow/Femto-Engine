#include "../f_ct.fh"


!  8 8888888888   8 8888888888            ,8.       ,8.    8888888 8888888888 ,o888888o.     
!  8 8888         8 8888                 ,888.     ,888.         8 8888    . 8888     `88.   
!  8 8888         8 8888                .`8888.   .`8888.        8 8888   ,8 8888       `8b  
!  8 8888         8 8888               ,8.`8888. ,8.`8888.       8 8888   88 8888        `8b 
!  8 888888888888 8 888888888888      ,8'8.`8888,8^8.`8888.      8 8888   88 8888         88 
!  8 8888         8 8888             ,8' `8.`8888' `8.`8888.     8 8888   88 8888         88 
!  8 8888         8 8888            ,8'   `8.`88'   `8.`8888.    8 8888   88 8888        ,8P 
!  8 8888         8 8888           ,8'     `8.`'     `8.`8888.   8 8888   `8 8888       ,8P  
!  8 8888         8 8888          ,8'       `8        `8.`8888.  8 8888    ` 8888     ,88'   
!  8 8888         8 888888888888 ,8'         `         `8.`8888. 8 8888       `8888888P'     



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y0(h, Y0, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: h(*), Y0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call g_diag_oovv_y0(h1, Y0, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_diag_oovv_y0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y0(h_, Y0_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y0_
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1
! Y0()  <-- 
! (    1.00000000)  h(c1,c1) 
do s_c1 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == 0) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Y0_ = Y0_ &
  + 1.00000000d+00 & 
  * h_(s_c1, s_c1)%array(i_c1, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x0(sc, ic, Y0, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y0
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x0(sc, ic, Y0, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x0(s_c, i_c, Y0, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y0
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (    2.00000000) Y0 D2(i,i,k,k) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * Y0 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x1(sc, ic, h, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: h(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x1(sc, ic, h1, av2_i2, d2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x1(s_c, i_c, h_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (    1.00000000) D2(i,i,k,k) h(a,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k) .and. &
IEOR(s_a,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i) & 
  * h_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x2(sc, ic, h, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: h(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x2(sc, ic, h1, av2_i2, d2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x2(s_c, i_c, h_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (    1.00000000) D2(i,i,k,k) h(c,c) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k) .and. &
IEOR(s_c,s_c) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i) & 
  * h_(s_c, s_c)%array(i_c, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y3(sc1, ic1, V2, Y1, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_diag_oovv_y3(sc1, ic1, h2_i, Y1, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_diag_oovv_y3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y3(s_c1, i_c1, V2_, Y1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y1_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y1()  <-- 
! (    1.00000000)  V2(c1,c1,c2,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == IEOR(s_c2,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y1_ = Y1_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c2, s_c1)%array(i_c2, i_c2, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x3(sc, ic, Y1, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y1
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x3(sc, ic, Y1, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x3(s_c, i_c, Y1, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y1
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (    2.00000000) Y1 D2(i,i,k,k) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * Y1 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y4(sc1, ic1, V2, Y2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y2
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_diag_oovv_y4(sc1, ic1, h2_i, Y2, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_diag_oovv_y4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y4(s_c1, i_c1, V2_, Y2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y2_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y2()  <-- 
! (    1.00000000)  V2(c1,c2,c1,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c2) == IEOR(s_c1,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y2_ = Y2_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c1, s_c2)%array(i_c2, i_c1, i_c2)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x4(sc, ic, Y2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y2
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x4(sc, ic, Y2, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x4(s_c, i_c, Y2, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y2
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (   -1.00000000) Y2 D2(i,i,k,k) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * Y2 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y5(sc1, ic1, V2, Y3, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y3, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_oovv_y5(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_oovv_y5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y5(s_c1, i_c1, V2_, Y3_, nir, nsym, psym)

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

integer :: s_a, i_a
! Y3(a,a)  <-- 
! (    1.00000000)  V2(c1,c1,a,a) 
do s_a = 0, nir-1
if( &
IEOR(s_a,s_a) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_a)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y3_(s_a, s_a)%array(i_a, i_a) =  &
    Y3_(s_a, s_a)%array(i_a, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_a, s_c1)%array(i_a, i_a, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x5(sc, ic, Y3, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y3(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y3, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x5(sc, ic, Yvv, av2_i2, d2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x5(s_c, i_c, Y3_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y3_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (    2.00000000) D2(i,i,k,k) Y3(a,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k) .and. &
IEOR(s_a,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i) & 
  * Y3_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y6(sc1, ic1, V2, Y4, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y4, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_oovv_y6(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_oovv_y6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y6(s_c1, i_c1, V2_, Y4_, nir, nsym, psym)

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

integer :: s_c, i_c
! Y4(c,c)  <-- 
! (    1.00000000)  V2(c1,c1,c,c) 
do s_c = 0, nir-1
if( &
IEOR(s_c,s_c) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c,s_c)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y4_(s_c, s_c)%array(i_c, i_c) =  &
    Y4_(s_c, s_c)%array(i_c, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_c, s_c1)%array(i_c, i_c, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x6(sc, ic, Y4, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y4(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y4, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x6(sc, ic, Yvv, av2_i2, d2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x6(s_c, i_c, Y4_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y4_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (    2.00000000) D2(i,i,k,k) Y4(c,c) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k) .and. &
IEOR(s_c,s_c) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i) & 
  * Y4_(s_c, s_c)%array(i_c, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y7(sc1, ic1, V2, Y5, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y5, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_oovv_y7(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_oovv_y7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y7(s_c1, i_c1, V2_, Y5_, nir, nsym, psym)

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

integer :: s_a, i_a
! Y5(a,a)  <-- 
! (    1.00000000)  V2(c1,a,c1,a) 
do s_a = 0, nir-1
if( &
IEOR(s_a,s_a) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_a)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y5_(s_a, s_a)%array(i_a, i_a) =  &
    Y5_(s_a, s_a)%array(i_a, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_c1, s_a)%array(i_a, i_c1, i_a)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x7(sc, ic, Y5, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y5(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y5, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x7(sc, ic, Yvv, av2_i2, d2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x7(s_c, i_c, Y5_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y5_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (   -1.00000000) D2(i,i,k,k) Y5(a,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k) .and. &
IEOR(s_a,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i) & 
  * Y5_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y8(sc1, ic1, V2, Y6, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y6, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_oovv_y8(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_oovv_y8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y8(s_c1, i_c1, V2_, Y6_, nir, nsym, psym)

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

integer :: s_c, i_c
! Y6(c,c)  <-- 
! (    1.00000000)  V2(c1,c,c1,c) 
do s_c = 0, nir-1
if( &
IEOR(s_c,s_c) == 0 .and. & 
IEOR(s_c1,s_c) == IEOR(s_c1,s_c)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y6_(s_c, s_c)%array(i_c, i_c) =  &
    Y6_(s_c, s_c)%array(i_c, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_c1, s_c)%array(i_c, i_c1, i_c)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x8(sc, ic, Y6, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y6(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y6, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x8(sc, ic, Yvv, av2_i2, d2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x8(s_c, i_c, Y6_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y6_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (   -1.00000000) D2(i,i,k,k) Y6(c,c) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k) .and. &
IEOR(s_c,s_c) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i) & 
  * Y6_(s_c, s_c)%array(i_c, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x9(sc, ic, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x9(sc, ic, h2_i, av2_i2, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x9(s_c, i_c, V2_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (    1.00000000) D2(i,i,k,k) V2(c,c,a,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k) .and. &
IEOR(s_c,s_c) == IEOR(s_a,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i) & 
  * V2_(s_a, s_a, s_c)%array(i_a, i_a, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y10(h, Y7, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: h(*), Y7
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call g_diag_oovv_y10(h1, Y7, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_diag_oovv_y10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y10(h_, Y7_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y7_
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1
! Y7()  <-- 
! (    1.00000000)  h(c1,c1) 
do s_c1 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == 0) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Y7_ = Y7_ &
  + 1.00000000d+00 & 
  * h_(s_c1, s_c1)%array(i_c1, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x10(sa, ia, Y7, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y7
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x10(sa, ia, Y7, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x10(s_a, i_a, Y7, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y7
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k
! Hdiag(i,k,a,a)  <-- 
! (    2.00000000) Y7 D2(i,k,k,i) 
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(s_i,s_k) == IEOR(s_k,s_i)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * Y7 & 
  * D2_(s_i, s_k, s_k, s_i)%array(i_i, i_k, i_k, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x11(sc, ic, h, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: h(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x11(sc, ic, h1, av2_i2, d3, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x11(s_c, i_c, h_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (    1.00000000) D3(i,i,k,k,o1,o2) h(o1,o2) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_i,s_i),s_k) == IEOR(IEOR(s_k,s_o1),s_o2) .and. &
IEOR(s_o1,s_o2) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o1, s_k, s_k, s_i, s_i)%array(i_o2, i_o1, i_k, i_k, i_i, i_i) & 
  * h_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_oovv_no0_x11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x12(sa, ia, h, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: h(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x12(sa, ia, h1, av2_i2, d3, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x12(s_a, i_a, h_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2
! Hdiag(i,k,a,a)  <-- 
! (    1.00000000) D3(i,k,k,i,o1,o2) h(o1,o2) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(IEOR(s_i,s_k),s_k) == IEOR(IEOR(s_i,s_o1),s_o2) .and. &
IEOR(s_o1,s_o2) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o1, s_i, s_k, s_k, s_i)%array(i_o2, i_o1, i_i, i_k, i_k, i_i) & 
  * h_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_oovv_no0_x12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x13(sa, ia, h, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: h(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x13(sa, ia, h1, av2_i2, d2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x13(s_a, i_a, h_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k
! Hdiag(i,k,a,a)  <-- 
! (    2.00000000) D2(i,k,k,i) h(a,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(s_i,s_k) == IEOR(s_k,s_i) .and. &
IEOR(s_a,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * D2_(s_i, s_k, s_k, s_i)%array(i_i, i_k, i_k, i_i) & 
  * h_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y14(sc1, ic1, V2, Y8, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y8
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_diag_oovv_y14(sc1, ic1, h2_i, Y8, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_diag_oovv_y14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y14(s_c1, i_c1, V2_, Y8_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y8_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y8()  <-- 
! (    1.00000000)  V2(c1,c1,c2,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == IEOR(s_c2,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y8_ = Y8_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c2, s_c1)%array(i_c2, i_c2, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x14(sa, ia, Y8, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y8
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x14(sa, ia, Y8, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x14(s_a, i_a, Y8, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y8
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k
! Hdiag(i,k,a,a)  <-- 
! (    2.00000000) Y8 D2(i,k,k,i) 
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(s_i,s_k) == IEOR(s_k,s_i)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * Y8 & 
  * D2_(s_i, s_k, s_k, s_i)%array(i_i, i_k, i_k, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y15(sc1, ic1, V2, Y9, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y9
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_diag_oovv_y15(sc1, ic1, h2_i, Y9, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_diag_oovv_y15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y15(s_c1, i_c1, V2_, Y9_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y9_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y9()  <-- 
! (    1.00000000)  V2(c1,c2,c1,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c2) == IEOR(s_c1,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y9_ = Y9_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c1, s_c2)%array(i_c2, i_c1, i_c2)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x15(sa, ia, Y9, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y9
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x15(sa, ia, Y9, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x15(s_a, i_a, Y9, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y9
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k
! Hdiag(i,k,a,a)  <-- 
! (   -1.00000000) Y9 D2(i,k,k,i) 
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(s_i,s_k) == IEOR(s_k,s_i)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * Y9 & 
  * D2_(s_i, s_k, s_k, s_i)%array(i_i, i_k, i_k, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y16(sc1, ic1, V2, Y10, nir, nsym, psym)

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
call g_diag_oovv_y16(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_oovv_y16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y16(s_c1, i_c1, V2_, Y10_, nir, nsym, psym)

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

end subroutine g_diag_oovv_y16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x16(sc, ic, Y10, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y10(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y10, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x16(sc, ic, Yaa, av2_i2, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x16(s_c, i_c, Y10_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y10_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (    2.00000000) D3(i,i,k,k,o1,o2) Y10(o1,o2) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_i,s_i),s_k) == IEOR(IEOR(s_k,s_o1),s_o2) .and. &
IEOR(s_o1,s_o2) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * D3_(s_o2, s_o1, s_k, s_k, s_i, s_i)%array(i_o2, i_o1, i_k, i_k, i_i, i_i) & 
  * Y10_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_oovv_no0_x16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y17(sc1, ic1, V2, Y11, nir, nsym, psym)

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
call g_diag_oovv_y17(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_oovv_y17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y17(s_c1, i_c1, V2_, Y11_, nir, nsym, psym)

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

end subroutine g_diag_oovv_y17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x17(sc, ic, Y11, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y11(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y11, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x17(sc, ic, Yaa, av2_i2, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x17(s_c, i_c, Y11_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y11_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (   -1.00000000) D3(i,i,k,k,o1,o2) Y11(o1,o2) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_i,s_i),s_k) == IEOR(IEOR(s_k,s_o1),s_o2) .and. &
IEOR(s_o1,s_o2) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_o2, s_o1, s_k, s_k, s_i, s_i)%array(i_o2, i_o1, i_k, i_k, i_i, i_i) & 
  * Y11_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_oovv_no0_x17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x18(sa, ia, sc, ic, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x18(sa, ia, sc, ic, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x18(s_a, i_a, s_c, i_c, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2
! Hdiag(i,k,a,c)  <-- 
! (    1.00000000) D3(i,i,k,k,o1,o2) V2(a,a,o1,o2) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_i,s_i),s_k) == IEOR(IEOR(s_k,s_o1),s_o2) .and. &
IEOR(s_a,s_a) == IEOR(s_o1,s_o2)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o1, s_k, s_k, s_i, s_i)%array(i_o2, i_o1, i_k, i_k, i_i, i_i) & 
  * V2_(s_o2, s_o1, s_a)%array(i_o2, i_o1, i_a)
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

end subroutine g_diag_oovv_no0_x18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x19(sc, ic, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x19(sc, ic, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x19(s_c, i_c, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (    1.00000000) D3(i,i,k,k,o1,o2) V2(c,c,o1,o2) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_i,s_i),s_k) == IEOR(IEOR(s_k,s_o1),s_o2) .and. &
IEOR(s_c,s_c) == IEOR(s_o1,s_o2)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o1, s_k, s_k, s_i, s_i)%array(i_o2, i_o1, i_k, i_k, i_i, i_i) & 
  * V2_(s_o2, s_o1, s_c)%array(i_o2, i_o1, i_c)
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

end subroutine g_diag_oovv_no0_x19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y20(sc1, ic1, V2, Y12, nir, nsym, psym)

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
call set_symblock_Yaa(sleft2, Y12, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_oovv_y20(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_oovv_y20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y20(s_c1, i_c1, V2_, Y12_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o2, i_o2
! Y12(o1,o2)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y12_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y12_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_o1, s_c1)%array(i_o2, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x20(sa, ia, Y12, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y12(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y12, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x20(sa, ia, Yaa, av2_i2, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x20(s_a, i_a, Y12_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y12_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2
! Hdiag(i,k,a,a)  <-- 
! (    1.00000000) D3(i,k,k,i,o1,o2) Y12(o1,o2) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(IEOR(s_i,s_k),s_k) == IEOR(IEOR(s_i,s_o1),s_o2) .and. &
IEOR(s_o1,s_o2) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o1, s_i, s_k, s_k, s_i)%array(i_o2, i_o1, i_i, i_k, i_k, i_i) & 
  * Y12_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_oovv_no0_x20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y21(sc1, ic1, V2, Y13, nir, nsym, psym)

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
call set_symblock_Yaa(sleft2, Y13, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_oovv_y21(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_oovv_y21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y21(s_c1, i_c1, V2_, Y13_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o2, i_o2
! Y13(o1,o2)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y13_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y13_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_c1, s_o1)%array(i_o2, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x21(sa, ia, Y13, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y13(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y13, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x21(sa, ia, Yaa, av2_i2, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x21(s_a, i_a, Y13_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y13_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2
! Hdiag(i,k,a,a)  <-- 
! (   -0.50000000) D3(i,k,k,i,o1,o2) Y13(o1,o2) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(IEOR(s_i,s_k),s_k) == IEOR(IEOR(s_i,s_o1),s_o2) .and. &
IEOR(s_o1,s_o2) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 0.50000000d+00 & 
  * D3_(s_o2, s_o1, s_i, s_k, s_k, s_i)%array(i_o2, i_o1, i_i, i_k, i_k, i_i) & 
  * Y13_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_oovv_no0_x21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y22(sc1, ic1, V2, Y14, nir, nsym, psym)

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
call g_diag_oovv_y22(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_oovv_y22



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y22(s_c1, i_c1, V2_, Y14_, nir, nsym, psym)

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

integer :: s_a, i_a, s_c, i_c
! Y14(a,c)  <-- 
! (    1.00000000)  V2(c1,c1,a,c) 
do s_a = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a,s_c) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_c)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y14_(s_c, s_a)%array(i_c, i_a) =  &
    Y14_(s_c, s_a)%array(i_c, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_a, s_c1)%array(i_c, i_a, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y22



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x22(sa, ia, Y14, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y14(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y14, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x22(sa, ia, Yvv, av2_i2, d2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x22



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x22(s_a, i_a, Y14_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y14_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k
! Hdiag(i,k,a,a)  <-- 
! (    2.00000000) D2(i,k,k,i) Y14(a,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(s_i,s_k) == IEOR(s_k,s_i) .and. &
IEOR(s_a,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * D2_(s_i, s_k, s_k, s_i)%array(i_i, i_k, i_k, i_i) & 
  * Y14_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x22



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y23(sc1, ic1, V2, Y15, nir, nsym, psym)

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
call g_diag_oovv_y23(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_oovv_y23



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y23(s_c1, i_c1, V2_, Y15_, nir, nsym, psym)

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

integer :: s_a, i_a, s_c, i_c
! Y15(a,c)  <-- 
! (    1.00000000)  V2(c1,a,c1,c) 
do s_a = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a,s_c) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_c)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y15_(s_c, s_a)%array(i_c, i_a) =  &
    Y15_(s_c, s_a)%array(i_c, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_c1, s_a)%array(i_c, i_c1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y23



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x23(sa, ia, Y15, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y15(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y15, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x23(sa, ia, Yvv, av2_i2, d2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x23



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x23(s_a, i_a, Y15_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y15_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k
! Hdiag(i,k,a,a)  <-- 
! (   -1.00000000) D2(i,k,k,i) Y15(a,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(s_i,s_k) == IEOR(s_k,s_i) .and. &
IEOR(s_a,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_i, s_k, s_k, s_i)%array(i_i, i_k, i_k, i_i) & 
  * Y15_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x23



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y24(sc1, ic1, V2, Y16, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y16(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y16, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_oovv_y24(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_oovv_y24



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y24(s_c1, i_c1, V2_, Y16_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y16_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y16(o1,o2)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y16_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y16_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_c1, s_o1)%array(i_o2, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y24



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x24(sa, ia, Y16, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y16(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y16, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x24(sa, ia, Yaa, av2_i2, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x24



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x24(s_a, i_a, Y16_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y16_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2
! Hdiag(i,k,a,a)  <-- 
! (   -0.50000000) D3(i,k,k,i,o1,o2) Y16(o1,o2) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(IEOR(s_i,s_k),s_k) == IEOR(IEOR(s_i,s_o1),s_o2) .and. &
IEOR(s_o1,s_o2) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 0.50000000d+00 & 
  * D3_(s_o2, s_o1, s_i, s_k, s_k, s_i)%array(i_o2, i_o1, i_i, i_k, i_k, i_i) & 
  * Y16_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_oovv_no0_x24



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y25(sc1, ic1, V2, Y17, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y17(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y17, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_oovv_y25(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_oovv_y25



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y25(s_c1, i_c1, V2_, Y17_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y17_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y17(o1,o2)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y17_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y17_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_o1, s_c1)%array(i_o2, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y25



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x25(sa, ia, Y17, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y17(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y17, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x25(sa, ia, Yaa, av2_i2, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x25



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x25(s_a, i_a, Y17_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y17_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2
! Hdiag(i,k,a,a)  <-- 
! (    1.00000000) D3(i,k,k,i,o1,o2) Y17(o1,o2) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(IEOR(s_i,s_k),s_k) == IEOR(IEOR(s_i,s_o1),s_o2) .and. &
IEOR(s_o1,s_o2) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o1, s_i, s_k, s_k, s_i)%array(i_o2, i_o1, i_i, i_k, i_k, i_i) & 
  * Y17_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_oovv_no0_x25



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x26(sc, ic, so1, io1, so3, io3, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, so1, io1, so3, io3
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x26(sc, ic, so1, io1, so3, io3, h2_i, av2_i2, d4_ij, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x26



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x26(s_c, i_c, s_o1, i_o1, s_o3, i_o3, V2_, Hdiag_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
integer, intent(in) :: i_o3, s_o3
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o2, i_o2, s_o4, i_o4, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (    0.50000000) D4(o1,o3,i,i,k,k,o2,o4) V2(o1,o3,o2,o4) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_o1,s_o3),IEOR(s_i,s_i)) == IEOR(IEOR(s_k,s_k),IEOR(s_o2,s_o4)) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 0.50000000d+00 & 
  * D4_(s_o4, s_o2, s_k, s_k, s_i, s_i)%array(i_o4, i_o2, i_k, i_k, i_i, i_i) & 
  * V2_(s_o4, s_o2, s_o3)%array(i_o4, i_o2, i_o3)
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

end subroutine g_diag_oovv_no0_x26



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x27(sa, ia, so1, io1, so3, io3, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, so1, io1, so3, io3
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x27(sa, ia, so1, io1, so3, io3, h2_i, av2_i2, d4_ij, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x27



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x27(s_a, i_a, s_o1, i_o1, s_o3, i_o3, V2_, Hdiag_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
integer, intent(in) :: i_o3, s_o3
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o2, i_o2, s_o4, i_o4
! Hdiag(i,k,a,a)  <-- 
! (    0.50000000) D4(o1,o3,i,k,k,i,o2,o4) V2(o1,o3,o2,o4) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(IEOR(s_o1,s_o3),IEOR(s_i,s_k)) == IEOR(IEOR(s_k,s_i),IEOR(s_o2,s_o4)) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 0.50000000d+00 & 
  * D4_(s_o4, s_o2, s_i, s_k, s_k, s_i)%array(i_o4, i_o2, i_i, i_k, i_k, i_i) & 
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

end subroutine g_diag_oovv_no0_x27



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x28(sa, ia, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x28(sa, ia, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x28



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x28(s_a, i_a, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2
! Hdiag(i,k,a,a)  <-- 
! (    1.00000000) D3(i,k,k,i,o1,o2) V2(a,a,o1,o2) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(IEOR(s_i,s_k),s_k) == IEOR(IEOR(s_i,s_o1),s_o2) .and. &
IEOR(s_a,s_a) == IEOR(s_o1,s_o2)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o1, s_i, s_k, s_k, s_i)%array(i_o2, i_o1, i_i, i_k, i_k, i_i) & 
  * V2_(s_o2, s_o1, s_a)%array(i_o2, i_o1, i_a)
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

end subroutine g_diag_oovv_no0_x28



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x29(sa, ia, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x29(sa, ia, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x29



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x29(s_a, i_a, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2
! Hdiag(i,k,a,a)  <-- 
! (    1.00000000) D3(i,k,k,o1,o2,i) V2(a,o1,o2,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(IEOR(s_i,s_k),s_k) == IEOR(IEOR(s_o1,s_o2),s_i) .and. &
IEOR(s_a,s_o1) == IEOR(s_o2,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_i, s_o2, s_o1, s_k, s_k, s_i)%array(i_i, i_o2, i_o1, i_k, i_k, i_i) & 
  * V2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
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

end subroutine g_diag_oovv_no0_x29



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x30(sa, ia, sc, ic, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x30(sa, ia, sc, ic, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x30



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x30(s_a, i_a, s_c, i_c, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2
! Hdiag(i,k,a,c)  <-- 
! (    1.00000000) D3(i,o1,k,k,o2,i) V2(a,o1,o2,a) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_i,s_o1),s_k) == IEOR(IEOR(s_k,s_o2),s_i) .and. &
IEOR(s_a,s_o1) == IEOR(s_o2,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_i, s_o2, s_k, s_k, s_o1, s_i)%array(i_i, i_o2, i_k, i_k, i_o1, i_i) & 
  * V2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
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

end subroutine g_diag_oovv_no0_x30



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x31(sc, ic, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x31(sc, ic, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x31



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x31(s_c, i_c, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (    1.00000000) D3(i,i,k,o1,o2,k) V2(c,o1,o2,c) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_i,s_i),s_k) == IEOR(IEOR(s_o1,s_o2),s_k) .and. &
IEOR(s_c,s_o1) == IEOR(s_o2,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_k, s_o2, s_o1, s_k, s_i, s_i)%array(i_k, i_o2, i_o1, i_k, i_i, i_i) & 
  * V2_(s_c, s_o2, s_o1)%array(i_c, i_o2, i_o1)
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

end subroutine g_diag_oovv_no0_x31



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y32(sc1, ic1, V2, Y18, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y18(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y18, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_oovv_y32(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_oovv_y32



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y32(s_c1, i_c1, V2_, Y18_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y18_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_c, i_c
! Y18(a,c)  <-- 
! (    1.00000000)  V2(c1,a,c1,c) 
do s_a = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a,s_c) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_c)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y18_(s_c, s_a)%array(i_c, i_a) =  &
    Y18_(s_c, s_a)%array(i_c, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_c1, s_a)%array(i_c, i_c1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y32



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x32(sa, ia, Y18, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y18(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y18, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x32(sa, ia, Yvv, av2_i2, d2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x32



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x32(s_a, i_a, Y18_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y18_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k
! Hdiag(i,k,a,a)  <-- 
! (   -1.00000000) D2(i,k,k,i) Y18(a,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(s_i,s_k) == IEOR(s_k,s_i) .and. &
IEOR(s_a,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_i, s_k, s_k, s_i)%array(i_i, i_k, i_k, i_i) & 
  * Y18_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x32



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_y33(sc1, ic1, V2, Y19, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y19(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y19, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_oovv_y33(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_oovv_y33



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_y33(s_c1, i_c1, V2_, Y19_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y19_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_c, i_c
! Y19(a,c)  <-- 
! (    1.00000000)  V2(c1,c1,a,c) 
do s_a = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a,s_c) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_c)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y19_(s_c, s_a)%array(i_c, i_a) =  &
    Y19_(s_c, s_a)%array(i_c, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_a, s_c1)%array(i_c, i_a, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_y33



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x33(sa, ia, Y19, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y19(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y19, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x33(sa, ia, Yvv, av2_i2, d2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x33



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x33(s_a, i_a, Y19_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y19_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k
! Hdiag(i,k,a,a)  <-- 
! (    2.00000000) D2(i,k,k,i) Y19(a,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(s_i,s_k) == IEOR(s_k,s_i) .and. &
IEOR(s_a,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * D2_(s_i, s_k, s_k, s_i)%array(i_i, i_k, i_k, i_i) & 
  * Y19_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x33



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x34(sa, ia, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x34(sa, ia, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x34



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x34(s_a, i_a, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2
! Hdiag(i,k,a,a)  <-- 
! (    1.00000000) D3(i,k,k,o1,o2,i) V2(a,o1,o2,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(IEOR(s_i,s_k),s_k) == IEOR(IEOR(s_o1,s_o2),s_i) .and. &
IEOR(s_a,s_o1) == IEOR(s_o2,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_i, s_o2, s_o1, s_k, s_k, s_i)%array(i_i, i_o2, i_o1, i_k, i_k, i_i) & 
  * V2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
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

end subroutine g_diag_oovv_no0_x34



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x35(sa, ia, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x35(sa, ia, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x35



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x35(s_a, i_a, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2
! Hdiag(i,k,a,a)  <-- 
! (    1.00000000) D3(i,k,k,i,o1,o2) V2(a,a,o1,o2) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_a) .and. & 
IEOR(IEOR(s_i,s_k),s_k) == IEOR(IEOR(s_i,s_o1),s_o2) .and. &
IEOR(s_a,s_a) == IEOR(s_o1,s_o2)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o1, s_i, s_k, s_k, s_i)%array(i_o2, i_o1, i_i, i_k, i_k, i_i) & 
  * V2_(s_o2, s_o1, s_a)%array(i_o2, i_o1, i_a)
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

end subroutine g_diag_oovv_no0_x35



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_no0_x36(sc, ic, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_no0_x36(sc, ic, h2_i, av2_i2, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_oovv_no0_x36



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_oovv_no0_x36(s_c, i_c, V2_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_a, i_a
! Hdiag(i,k,a,c)  <-- 
! (    1.00000000) D2(i,k,k,i) V2(c,a,a,c) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_k) == IEOR(s_k,s_i) .and. &
IEOR(s_c,s_a) == IEOR(s_a,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    Hdiag_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_i, s_k, s_k, s_i)%array(i_i, i_k, i_k, i_i) & 
  * V2_(s_c, s_a, s_a)%array(i_c, i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_oovv_no0_x36

