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
subroutine g_if_diag_ccvv_y0(h, Y0, nir, nsym, psym)

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
call g_diag_ccvv_y0(h1, Y0, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_diag_ccvv_y0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y0(h_, Y0_, nir, nsym, psym)

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

end subroutine g_diag_ccvv_y0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x0(sc, ic, Y0, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x0(sc, ic, Y0, av2_i2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x0(s_c, i_c, Y0, Hdiag_, nir, nsym, psym)

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

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    8.00000000) Y0 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * Y0
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y1(h, Y1, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: h(*), Y1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call g_diag_ccvv_y1(h1, Y1, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_diag_ccvv_y1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y1(h_, Y1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y1_
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1
! Y1()  <-- 
! (    1.00000000)  h(c1,c1) 
do s_c1 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == 0) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Y1_ = Y1_ &
  + 1.00000000d+00 & 
  * h_(s_c1, s_c1)%array(i_c1, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x1(sc, ic, Y1, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x1(sc, ic, Y1, av2_i2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x1(s_c, i_c, Y1, Hdiag_, nir, nsym, psym)

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

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (   -4.00000000) Y1 
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 4.00000000d+00 & 
  * Y1
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x2(sc, ic, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x2(sc, ic, h1, av2_i2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x2(s_c, i_c, h_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (   -4.00000000) h(w,w) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_w) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * h_(s_w, s_w)%array(i_w, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x3(sc, ic, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x3(sc, ic, h1, av2_i2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x3(s_c, i_c, h_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (    4.00000000) h(w,w) 
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_w) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 4.00000000d+00 & 
  * h_(s_w, s_w)%array(i_w, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x4(sc, ic, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x4(sc, ic, h1, av2_i2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x4(s_c, i_c, h_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (   -4.00000000) h(y,y) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_y) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * h_(s_y, s_y)%array(i_y, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y5(h, Y2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: h(*), Y2
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call g_diag_ccvv_y5(h1, Y2, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_diag_ccvv_y5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y5(h_, Y2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y2_
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1
! Y2()  <-- 
! (    1.00000000)  h(c1,c1) 
do s_c1 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == 0) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Y2_ = Y2_ &
  + 1.00000000d+00 & 
  * h_(s_c1, s_c1)%array(i_c1, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x5(sa, ia, Y2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y2
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x5(sa, ia, Y2, av2_i2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x5(s_a, i_a, Y2, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y2
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    8.00000000) Y2 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 8.00000000d+00 & 
  * Y2
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y6(h, Y3, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: h(*), Y3
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call g_diag_ccvv_y6(h1, Y3, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_diag_ccvv_y6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y6(h_, Y3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y3_
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1
! Y3()  <-- 
! (    1.00000000)  h(c1,c1) 
do s_c1 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == 0) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Y3_ = Y3_ &
  + 1.00000000d+00 & 
  * h_(s_c1, s_c1)%array(i_c1, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x6(sa, ia, Y3, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y3
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x6(sa, ia, Y3, av2_i2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x6(s_a, i_a, Y3, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y3
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -4.00000000) Y3 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * Y3
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x7(sa, ia, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x7(sa, ia, h1, av2_i2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x7(s_a, i_a, h_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (   -8.00000000) h(w,w) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_w,s_w) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 8.00000000d+00 & 
  * h_(s_w, s_w)%array(i_w, i_w)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x8(sa, ia, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x8(sa, ia, h1, av2_i2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x8(s_a, i_a, h_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (    2.00000000) h(w,w) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_w,s_w) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * h_(s_w, s_w)%array(i_w, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x9(sa, ia, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x9(sa, ia, h1, av2_i2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x9(s_a, i_a, h_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! Hdiag(w,y,a,a)  <-- 
! (    2.00000000) h(y,y) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_y,s_y) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * h_(s_y, s_y)%array(i_y, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x10(sc, ic, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x10(sc, ic, h1, av2_i2, d1, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x10(s_c, i_c, h_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    4.00000000) D1(o1,o2) h(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
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

end subroutine g_diag_ccvv_no0_x10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x11(sc, ic, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x11(sc, ic, h1, av2_i2, d1, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x11(s_c, i_c, h_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (   -2.00000000) D1(o1,o2) h(o2,o1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o2,s_o1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
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

end subroutine g_diag_ccvv_no0_x11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x12(sa, ia, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x12(sa, ia, h1, av2_i2, d1, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x12(s_a, i_a, h_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    4.00000000) D1(o1,o2) h(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * h_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x13(sa, ia, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x13(sa, ia, h1, av2_i2, d1, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x13(s_a, i_a, h_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -2.00000000) D1(o1,o2) h(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
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

end subroutine g_diag_ccvv_no0_x13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x14(sa, ia, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x14(sa, ia, h1, av2_i2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x14(s_a, i_a, h_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    8.00000000) h(a,a) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 8.00000000d+00 & 
  * h_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x15(sa, ia, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x15(sa, ia, h1, av2_i2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x15(s_a, i_a, h_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -4.00000000) h(a,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * h_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x16(sc, ic, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x16(sc, ic, h1, av2_i2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x16(s_c, i_c, h_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,c)  <-- 
! (    4.00000000) h(a,a) 
do s_a = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_a,s_a) == 0) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * h_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x17(sc, ic, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x17(sc, ic, h1, av2_i2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x17(s_c, i_c, h_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_w, i_w
! Hdiag(w,w,a,c)  <-- 
! (   -2.00000000) h(a,a) 
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_a,s_a) == 0) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 2.00000000d+00 & 
  * h_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x18(sc, ic, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x18(sc, ic, h1, av2_i2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x18(s_c, i_c, h_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (   -2.00000000) h(c,c) 
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_c) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 2.00000000d+00 & 
  * h_(s_c, s_c)%array(i_c, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y19(sc1, ic1, V2, Y4, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y4
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_diag_ccvv_y19(sc1, ic1, h2_i, Y4, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_diag_ccvv_y19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y19(s_c1, i_c1, V2_, Y4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y4_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y4()  <-- 
! (    1.00000000)  V2(c1,c1,c2,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == IEOR(s_c2,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y4_ = Y4_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c2, s_c1)%array(i_c2, i_c2, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x19(sc, ic, Y4, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y4
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x19(sc, ic, Y4, av2_i2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x19(s_c, i_c, Y4, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y4
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (   -4.00000000) Y4 
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 4.00000000d+00 & 
  * Y4
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y20(sc1, ic1, V2, Y5, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y5
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_diag_ccvv_y20(sc1, ic1, h2_i, Y5, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_diag_ccvv_y20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y20(s_c1, i_c1, V2_, Y5_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y5_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y5()  <-- 
! (    1.00000000)  V2(c1,c2,c1,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c2) == IEOR(s_c1,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y5_ = Y5_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c1, s_c2)%array(i_c2, i_c1, i_c2)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x20(sc, ic, Y5, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y5
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x20(sc, ic, Y5, av2_i2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x20(s_c, i_c, Y5, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y5
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (    2.00000000) Y5 
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 2.00000000d+00 & 
  * Y5
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y21(sc1, ic1, V2, Y6, nir, nsym, psym)

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
call g_diag_ccvv_y21(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y21(s_c1, i_c1, V2_, Y6_, nir, nsym, psym)

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

integer :: s_a, i_a
! Y6(a,a)  <-- 
! (    1.00000000)  V2(c1,c1,a,a) 
do s_a = 0, nir-1
if( &
IEOR(s_a,s_a) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_a)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y6_(s_a, s_a)%array(i_a, i_a) =  &
    Y6_(s_a, s_a)%array(i_a, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_a, s_c1)%array(i_a, i_a, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x21(sc, ic, Y6, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x21(sc, ic, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x21(s_c, i_c, Y6_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y6_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_w, i_w
! Hdiag(w,w,a,c)  <-- 
! (   -4.00000000) Y6(a,a) 
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_a,s_a) == 0) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 4.00000000d+00 & 
  * Y6_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y22(sc1, ic1, V2, Y7, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y7, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ccvv_y22(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y22



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y22(s_c1, i_c1, V2_, Y7_, nir, nsym, psym)

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

integer :: s_c, i_c
! Y7(c,c)  <-- 
! (    1.00000000)  V2(c1,c1,c,c) 
do s_c = 0, nir-1
if( &
IEOR(s_c,s_c) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c,s_c)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y7_(s_c, s_c)%array(i_c, i_c) =  &
    Y7_(s_c, s_c)%array(i_c, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_c, s_c1)%array(i_c, i_c, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y22



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x22(sc, ic, Y7, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y7(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y7, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x22(sc, ic, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x22



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x22(s_c, i_c, Y7_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y7_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (   -4.00000000) Y7(c,c) 
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_c) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 4.00000000d+00 & 
  * Y7_(s_c, s_c)%array(i_c, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x22



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y23(sc1, ic1, V2, Y8, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y8, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ccvv_y23(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y23



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y23(s_c1, i_c1, V2_, Y8_, nir, nsym, psym)

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

integer :: s_a, i_a
! Y8(a,a)  <-- 
! (    1.00000000)  V2(c1,a,c1,a) 
do s_a = 0, nir-1
if( &
IEOR(s_a,s_a) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_a)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y8_(s_a, s_a)%array(i_a, i_a) =  &
    Y8_(s_a, s_a)%array(i_a, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_c1, s_a)%array(i_a, i_c1, i_a)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y23



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x23(sc, ic, Y8, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y8(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y8, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x23(sc, ic, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x23



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x23(s_c, i_c, Y8_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y8_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_w, i_w
! Hdiag(w,w,a,c)  <-- 
! (    2.00000000) Y8(a,a) 
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_a,s_a) == 0) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 2.00000000d+00 & 
  * Y8_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x23



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y24(sc1, ic1, V2, Y9, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y9, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ccvv_y24(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y24



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y24(s_c1, i_c1, V2_, Y9_, nir, nsym, psym)

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

integer :: s_c, i_c
! Y9(c,c)  <-- 
! (    1.00000000)  V2(c1,c,c1,c) 
do s_c = 0, nir-1
if( &
IEOR(s_c,s_c) == 0 .and. & 
IEOR(s_c1,s_c) == IEOR(s_c1,s_c)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y9_(s_c, s_c)%array(i_c, i_c) =  &
    Y9_(s_c, s_c)%array(i_c, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_c1, s_c)%array(i_c, i_c1, i_c)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y24



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x24(sc, ic, Y9, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y9(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y9, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x24(sc, ic, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x24



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x24(s_c, i_c, Y9_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y9_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (    2.00000000) Y9(c,c) 
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_c) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 2.00000000d+00 & 
  * Y9_(s_c, s_c)%array(i_c, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x24



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x25(sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x25(sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x25



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x25(s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_w, i_w
! Hdiag(w,w,a,c)  <-- 
! (   -2.00000000) V2(c,c,a,a) 
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_c) == IEOR(s_a,s_a)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 2.00000000d+00 & 
  * V2_(s_a, s_a, s_c)%array(i_a, i_a, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x25



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x26(sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x26(sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x26



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x26(s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_w, i_w
! Hdiag(w,w,a,c)  <-- 
! (    4.00000000) V2(c,a,a,c) 
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_a) == IEOR(s_a,s_c)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 4.00000000d+00 & 
  * V2_(s_c, s_a, s_a)%array(i_c, i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x26



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x27(sc, ic, h, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x27(sc, ic, h1, av2_i2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x27



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x27(s_c, i_c, h_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    4.00000000) h(c,c) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_c) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * h_(s_c, s_c)%array(i_c, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x27



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y28(sc1, ic1, V2, Y10, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y10
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_diag_ccvv_y28(sc1, ic1, h2_i, Y10, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_diag_ccvv_y28



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y28(s_c1, i_c1, V2_, Y10_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y10_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y10()  <-- 
! (    1.00000000)  V2(c1,c1,c2,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == IEOR(s_c2,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y10_ = Y10_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c2, s_c1)%array(i_c2, i_c2, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y28



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x28(sc, ic, Y10, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y10
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x28(sc, ic, Y10, av2_i2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x28



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x28(s_c, i_c, Y10, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y10
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    8.00000000) Y10 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * Y10
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x28



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y29(sc1, ic1, V2, Y11, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y11, nir, nsym, psym) ! -> Ycc (allocate) 
call g_diag_ccvv_y29(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_diag_ccvv_y29



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y29(s_c1, i_c1, V2_, Y11_, nir, nsym, psym)

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

integer :: s_w, i_w
! Y11(w,w)  <-- 
! (    1.00000000)  V2(c1,c1,w,w) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_w,s_w)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Y11_(s_w, s_w)%array(i_w, i_w) =  &
    Y11_(s_w, s_w)%array(i_w, i_w) &
  + 1.00000000d+00 & 
  * V2_(s_w, s_w, s_c1)%array(i_w, i_w, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y29



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x29(sc, ic, Y11, Hdiag, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y11, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x29(sc, ic, Ycc, av2_i2, nir, nsym, psym)

deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x29



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x29(s_c, i_c, Y11_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y11_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (   -8.00000000) Y11(w,w) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_w) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 8.00000000d+00 & 
  * Y11_(s_w, s_w)%array(i_w, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x29



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y30(sc1, ic1, V2, Y12, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y12, nir, nsym, psym) ! -> Ycc (allocate) 
call g_diag_ccvv_y30(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_diag_ccvv_y30



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y30(s_c1, i_c1, V2_, Y12_, nir, nsym, psym)

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

integer :: s_w, i_w, s_y, i_y
! Y12(w,y)  <-- 
! (    1.00000000)  V2(c1,c1,w,y) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_w,s_y)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y12_(s_y, s_w)%array(i_y, i_w) =  &
    Y12_(s_y, s_w)%array(i_y, i_w) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_w, s_c1)%array(i_y, i_w, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y30



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x30(sc, ic, Y12, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y12(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Ycc(sleft2, Y12, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x30(sc, ic, Ycc, av2_i2, nir, nsym, psym)

deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x30



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x30(s_c, i_c, Y12_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y12_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (    8.00000000) Y12(w,w) 
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_w) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 8.00000000d+00 & 
  * Y12_(s_w, s_w)%array(i_w, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x30



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y31(sc1, ic1, V2, Y13, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y13, nir, nsym, psym) ! -> Ycc (allocate) 
call g_diag_ccvv_y31(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_diag_ccvv_y31



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y31(s_c1, i_c1, V2_, Y13_, nir, nsym, psym)

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

integer :: s_y, i_y
! Y13(y,y)  <-- 
! (    1.00000000)  V2(c1,c1,y,y) 
do s_y = 0, nir-1
if( &
IEOR(s_y,s_y) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_y,s_y)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y13_(s_y, s_y)%array(i_y, i_y) =  &
    Y13_(s_y, s_y)%array(i_y, i_y) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_y, s_c1)%array(i_y, i_y, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y31



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x31(sc, ic, Y13, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y13(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Ycc(sleft2, Y13, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x31(sc, ic, Ycc, av2_i2, nir, nsym, psym)

deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x31



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x31(s_c, i_c, Y13_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y13_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (   -8.00000000) Y13(y,y) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_y) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 8.00000000d+00 & 
  * Y13_(s_y, s_y)%array(i_y, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x31



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y32(sc1, ic1, V2, Y14, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y14
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_diag_ccvv_y32(sc1, ic1, h2_i, Y14, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_diag_ccvv_y32



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y32(s_c1, i_c1, V2_, Y14_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y14_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y14()  <-- 
! (    1.00000000)  V2(c1,c2,c1,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c2) == IEOR(s_c1,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y14_ = Y14_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c1, s_c2)%array(i_c2, i_c1, i_c2)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y32



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x32(sc, ic, Y14, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y14
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x32(sc, ic, Y14, av2_i2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x32



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x32(s_c, i_c, Y14, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y14
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (   -4.00000000) Y14 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * Y14
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x32



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y33(sc1, ic1, V2, Y15, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y15, nir, nsym, psym) ! -> Ycc (allocate) 
call g_diag_ccvv_y33(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_diag_ccvv_y33



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y33(s_c1, i_c1, V2_, Y15_, nir, nsym, psym)

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

integer :: s_w, i_w
! Y15(w,w)  <-- 
! (    1.00000000)  V2(c1,w,c1,w) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == 0 .and. & 
IEOR(s_c1,s_w) == IEOR(s_c1,s_w)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Y15_(s_w, s_w)%array(i_w, i_w) =  &
    Y15_(s_w, s_w)%array(i_w, i_w) &
  + 1.00000000d+00 & 
  * V2_(s_w, s_c1, s_w)%array(i_w, i_c1, i_w)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y33



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x33(sc, ic, Y15, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y15(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Ycc(sleft2, Y15, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x33(sc, ic, Ycc, av2_i2, nir, nsym, psym)

deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x33



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x33(s_c, i_c, Y15_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y15_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    4.00000000) Y15(w,w) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_w) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * Y15_(s_w, s_w)%array(i_w, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x33



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y34(sc1, ic1, V2, Y16, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y16, nir, nsym, psym) ! -> Ycc (allocate) 
call g_diag_ccvv_y34(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_diag_ccvv_y34



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y34(s_c1, i_c1, V2_, Y16_, nir, nsym, psym)

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

integer :: s_w, i_w, s_y, i_y
! Y16(w,y)  <-- 
! (    1.00000000)  V2(c1,w,c1,y) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == 0 .and. & 
IEOR(s_c1,s_w) == IEOR(s_c1,s_y)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y16_(s_y, s_w)%array(i_y, i_w) =  &
    Y16_(s_y, s_w)%array(i_y, i_w) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_c1, s_w)%array(i_y, i_c1, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y34



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x34(sc, ic, Y16, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y16(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Ycc(sleft2, Y16, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x34(sc, ic, Ycc, av2_i2, nir, nsym, psym)

deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x34



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x34(s_c, i_c, Y16_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y16_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (   -4.00000000) Y16(w,w) 
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_w) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 4.00000000d+00 & 
  * Y16_(s_w, s_w)%array(i_w, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x34



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y35(sc1, ic1, V2, Y17, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y17, nir, nsym, psym) ! -> Ycc (allocate) 
call g_diag_ccvv_y35(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_diag_ccvv_y35



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y35(s_c1, i_c1, V2_, Y17_, nir, nsym, psym)

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

integer :: s_y, i_y
! Y17(y,y)  <-- 
! (    1.00000000)  V2(c1,y,c1,y) 
do s_y = 0, nir-1
if( &
IEOR(s_y,s_y) == 0 .and. & 
IEOR(s_c1,s_y) == IEOR(s_c1,s_y)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y17_(s_y, s_y)%array(i_y, i_y) =  &
    Y17_(s_y, s_y)%array(i_y, i_y) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_c1, s_y)%array(i_y, i_c1, i_y)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y35



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x35(sc, ic, Y17, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y17(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Ycc(sleft2, Y17, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x35(sc, ic, Ycc, av2_i2, nir, nsym, psym)

deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x35



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x35(s_c, i_c, Y17_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y17_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    4.00000000) Y17(y,y) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_y) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * Y17_(s_y, s_y)%array(i_y, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x35



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x36(sc, ic, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x36(sc, ic, sw, iw, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x36



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x36(s_c, i_c, s_w, i_w, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    4.00000000) V2(w,w,y,y) 
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_w) == IEOR(s_y,s_y)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * V2_(s_y, s_y, s_w)%array(i_y, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x36



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x37(sc, ic, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x37(sc, ic, sw, iw, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x37



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x37(s_c, i_c, s_w, i_w, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (   -2.00000000) V2(w,y,w,y) 
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_w,s_y)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * V2_(s_y, s_w, s_y)%array(i_y, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x37



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y38(sc1, ic1, V2, Y18, nir, nsym, psym)

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
call g_diag_ccvv_y38(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y38



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y38(s_c1, i_c1, V2_, Y18_, nir, nsym, psym)

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

integer :: s_a, i_a
! Y18(a,a)  <-- 
! (    1.00000000)  V2(c1,c1,a,a) 
do s_a = 0, nir-1
if( &
IEOR(s_a,s_a) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_a)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y18_(s_a, s_a)%array(i_a, i_a) =  &
    Y18_(s_a, s_a)%array(i_a, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_a, s_c1)%array(i_a, i_a, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y38



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x38(sc, ic, Y18, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y18(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y18, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x38(sc, ic, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x38



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x38(s_c, i_c, Y18_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y18_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,c)  <-- 
! (    8.00000000) Y18(a,a) 
do s_a = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_a,s_a) == 0) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * Y18_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x38



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x39(sa, ia, sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x39(sa, ia, sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x39



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x39(s_a, i_a, s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,c)  <-- 
! (   -4.00000000) V2(a,a,w,w) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_a,s_a) == IEOR(s_w,s_w)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * V2_(s_w, s_w, s_a)%array(i_w, i_w, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x39



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x40(sa, ia, sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x40(sa, ia, sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x40



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x40(s_a, i_a, s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! Hdiag(w,y,a,c)  <-- 
! (   -4.00000000) V2(a,a,y,y) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_a,s_a) == IEOR(s_y,s_y)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * V2_(s_y, s_y, s_a)%array(i_y, i_y, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x40



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y41(sc1, ic1, V2, Y19, nir, nsym, psym)

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
call g_diag_ccvv_y41(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y41



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y41(s_c1, i_c1, V2_, Y19_, nir, nsym, psym)

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

integer :: s_c, i_c
! Y19(c,c)  <-- 
! (    1.00000000)  V2(c1,c1,c,c) 
do s_c = 0, nir-1
if( &
IEOR(s_c,s_c) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c,s_c)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y19_(s_c, s_c)%array(i_c, i_c) =  &
    Y19_(s_c, s_c)%array(i_c, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_c, s_c1)%array(i_c, i_c, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y41



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x41(sc, ic, Y19, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y19(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y19, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x41(sc, ic, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x41



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x41(s_c, i_c, Y19_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y19_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    8.00000000) Y19(c,c) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_c) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * Y19_(s_c, s_c)%array(i_c, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x41



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x42(sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x42(sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x42



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x42(s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (   -4.00000000) V2(c,c,w,w) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_c) == IEOR(s_w,s_w)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * V2_(s_w, s_w, s_c)%array(i_w, i_w, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x42



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x43(sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x43(sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x43



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x43(s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (   -4.00000000) V2(c,c,y,y) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_c) == IEOR(s_y,s_y)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * V2_(s_y, s_y, s_c)%array(i_y, i_y, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x43



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x44(sa, ia, sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x44(sa, ia, sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x44



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x44(s_a, i_a, s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,c)  <-- 
! (    8.00000000) V2(a,w,w,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_a,s_w) == IEOR(s_w,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * V2_(s_a, s_w, s_w)%array(i_a, i_w, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x44



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y45(sc1, ic1, V2, Y20, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y20(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y20, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ccvv_y45(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y45



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y45(s_c1, i_c1, V2_, Y20_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y20_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a
! Y20(a,a)  <-- 
! (    1.00000000)  V2(c1,a,c1,a) 
do s_a = 0, nir-1
if( &
IEOR(s_a,s_a) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_a)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y20_(s_a, s_a)%array(i_a, i_a) =  &
    Y20_(s_a, s_a)%array(i_a, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_c1, s_a)%array(i_a, i_c1, i_a)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y45



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x45(sc, ic, Y20, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y20(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y20, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x45(sc, ic, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x45



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x45(s_c, i_c, Y20_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y20_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,c)  <-- 
! (   -4.00000000) Y20(a,a) 
do s_a = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_a,s_a) == 0) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * Y20_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x45



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x46(sa, ia, sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x46(sa, ia, sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x46



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x46(s_a, i_a, s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! Hdiag(w,y,a,c)  <-- 
! (    2.00000000) V2(a,y,y,a) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_a,s_y) == IEOR(s_y,s_a)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * V2_(s_a, s_y, s_y)%array(i_a, i_y, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x46



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x47(sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x47(sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x47



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x47(s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    8.00000000) V2(c,y,y,c) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_y) == IEOR(s_y,s_c)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * V2_(s_c, s_y, s_y)%array(i_c, i_y, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x47



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x48(sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x48(sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x48



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x48(s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    2.00000000) V2(c,w,w,c) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_w) == IEOR(s_w,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * V2_(s_c, s_w, s_w)%array(i_c, i_w, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x48



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y49(sc1, ic1, V2, Y21, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y21(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y21, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ccvv_y49(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y49



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y49(s_c1, i_c1, V2_, Y21_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y21_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c, i_c
! Y21(c,c)  <-- 
! (    1.00000000)  V2(c1,c,c1,c) 
do s_c = 0, nir-1
if( &
IEOR(s_c,s_c) == 0 .and. & 
IEOR(s_c1,s_c) == IEOR(s_c1,s_c)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y21_(s_c, s_c)%array(i_c, i_c) =  &
    Y21_(s_c, s_c)%array(i_c, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_c1, s_c)%array(i_c, i_c1, i_c)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y49



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x49(sc, ic, Y21, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y21(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y21, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x49(sc, ic, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x49



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x49(s_c, i_c, Y21_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y21_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (   -4.00000000) Y21(c,c) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_c) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * Y21_(s_c, s_c)%array(i_c, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x49



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x50(sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x50(sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x50



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x50(s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,c)  <-- 
! (    4.00000000) V2(c,c,a,a) 
do s_a = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_c) == IEOR(s_a,s_a)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * V2_(s_a, s_a, s_c)%array(i_a, i_a, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x50



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x51(sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x51(sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x51



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x51(s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,c)  <-- 
! (   -2.00000000) V2(c,a,a,c) 
do s_a = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_a) == IEOR(s_a,s_c)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * V2_(s_c, s_a, s_a)%array(i_c, i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x51



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y52(sc1, ic1, V2, Y22, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y22
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_diag_ccvv_y52(sc1, ic1, h2_i, Y22, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_diag_ccvv_y52



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y52(s_c1, i_c1, V2_, Y22_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y22_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y22()  <-- 
! (    1.00000000)  V2(c1,c1,c2,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == IEOR(s_c2,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y22_ = Y22_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c2, s_c1)%array(i_c2, i_c2, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y52



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x52(sa, ia, Y22, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y22
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x52(sa, ia, Y22, av2_i2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x52



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x52(s_a, i_a, Y22, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y22
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    8.00000000) Y22 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 8.00000000d+00 & 
  * Y22
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x52



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y53(sc1, ic1, V2, Y23, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y23
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_diag_ccvv_y53(sc1, ic1, h2_i, Y23, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_diag_ccvv_y53



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y53(s_c1, i_c1, V2_, Y23_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y23_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y23()  <-- 
! (    1.00000000)  V2(c1,c1,c2,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == IEOR(s_c2,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y23_ = Y23_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c2, s_c1)%array(i_c2, i_c2, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y53



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x53(sa, ia, Y23, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y23
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x53(sa, ia, Y23, av2_i2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x53



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x53(s_a, i_a, Y23, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y23
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -4.00000000) Y23 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * Y23
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x53



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y54(sc1, ic1, V2, Y24, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y24(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y24, nir, nsym, psym) ! -> Ycc (allocate) 
call g_diag_ccvv_y54(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_diag_ccvv_y54



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y54(s_c1, i_c1, V2_, Y24_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y24_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Y24(w,y)  <-- 
! (    1.00000000)  V2(c1,c1,w,y) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_w,s_y)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y24_(s_y, s_w)%array(i_y, i_w) =  &
    Y24_(s_y, s_w)%array(i_y, i_w) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_w, s_c1)%array(i_y, i_w, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y54



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x54(sa, ia, Y24, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y24(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Ycc(sleft2, Y24, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x54(sa, ia, Ycc, av2_i2, nir, nsym, psym)

deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x54



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x54(s_a, i_a, Y24_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y24_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (  -16.00000000) Y24(w,w) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_w,s_w) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 16.00000000d+00 & 
  * Y24_(s_w, s_w)%array(i_w, i_w)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x54



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y55(sc1, ic1, V2, Y25, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y25(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y25, nir, nsym, psym) ! -> Ycc (allocate) 
call g_diag_ccvv_y55(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_diag_ccvv_y55



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y55(s_c1, i_c1, V2_, Y25_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y25_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Y25(w,w)  <-- 
! (    1.00000000)  V2(c1,c1,w,w) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_w,s_w)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Y25_(s_w, s_w)%array(i_w, i_w) =  &
    Y25_(s_w, s_w)%array(i_w, i_w) &
  + 1.00000000d+00 & 
  * V2_(s_w, s_w, s_c1)%array(i_w, i_w, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y55



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x55(sa, ia, Y25, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y25(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Ycc(sleft2, Y25, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x55(sa, ia, Ycc, av2_i2, nir, nsym, psym)

deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x55



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x55(s_a, i_a, Y25_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y25_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (    4.00000000) Y25(w,w) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_w,s_w) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * Y25_(s_w, s_w)%array(i_w, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x55



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y56(sc1, ic1, V2, Y26, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y26(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y26, nir, nsym, psym) ! -> Ycc (allocate) 
call g_diag_ccvv_y56(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_diag_ccvv_y56



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y56(s_c1, i_c1, V2_, Y26_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y26_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y
! Y26(y,y)  <-- 
! (    1.00000000)  V2(c1,c1,y,y) 
do s_y = 0, nir-1
if( &
IEOR(s_y,s_y) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_y,s_y)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y26_(s_y, s_y)%array(i_y, i_y) =  &
    Y26_(s_y, s_y)%array(i_y, i_y) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_y, s_c1)%array(i_y, i_y, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y56



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x56(sa, ia, Y26, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y26(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Ycc(sleft2, Y26, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x56(sa, ia, Ycc, av2_i2, nir, nsym, psym)

deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x56



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x56(s_a, i_a, Y26_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y26_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! Hdiag(w,y,a,a)  <-- 
! (    4.00000000) Y26(y,y) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_y,s_y) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * Y26_(s_y, s_y)%array(i_y, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x56



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y57(sc1, ic1, V2, Y27, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y27
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_diag_ccvv_y57(sc1, ic1, h2_i, Y27, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_diag_ccvv_y57



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y57(s_c1, i_c1, V2_, Y27_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y27_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y27()  <-- 
! (    1.00000000)  V2(c1,c2,c1,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c2) == IEOR(s_c1,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y27_ = Y27_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c1, s_c2)%array(i_c2, i_c1, i_c2)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y57



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x57(sa, ia, Y27, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y27
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x57(sa, ia, Y27, av2_i2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x57



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x57(s_a, i_a, Y27, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y27
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (   -4.00000000) Y27 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 4.00000000d+00 & 
  * Y27
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x57



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y58(sc1, ic1, V2, Y28, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y28
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_diag_ccvv_y58(sc1, ic1, h2_i, Y28, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_diag_ccvv_y58



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y58(s_c1, i_c1, V2_, Y28_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y28_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y28()  <-- 
! (    1.00000000)  V2(c1,c2,c1,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c2) == IEOR(s_c1,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y28_ = Y28_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c1, s_c2)%array(i_c2, i_c1, i_c2)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y58



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x58(sa, ia, Y28, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y28
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x58(sa, ia, Y28, av2_i2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x58



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x58(s_a, i_a, Y28, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y28
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (    2.00000000) Y28 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * Y28
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x58



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y59(sc1, ic1, V2, Y29, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y29(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y29, nir, nsym, psym) ! -> Ycc (allocate) 
call g_diag_ccvv_y59(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_diag_ccvv_y59



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y59(s_c1, i_c1, V2_, Y29_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y29_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Y29(w,w)  <-- 
! (    1.00000000)  V2(c1,w,c1,w) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == 0 .and. & 
IEOR(s_c1,s_w) == IEOR(s_c1,s_w)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Y29_(s_w, s_w)%array(i_w, i_w) =  &
    Y29_(s_w, s_w)%array(i_w, i_w) &
  + 1.00000000d+00 & 
  * V2_(s_w, s_c1, s_w)%array(i_w, i_c1, i_w)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y59



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x59(sa, ia, Y29, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y29(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Ycc(sleft2, Y29, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x59(sa, ia, Ycc, av2_i2, nir, nsym, psym)

deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x59



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x59(s_a, i_a, Y29_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y29_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -2.00000000) Y29(w,w) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_w,s_w) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * Y29_(s_w, s_w)%array(i_w, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x59



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y60(sc1, ic1, V2, Y30, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y30(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y30, nir, nsym, psym) ! -> Ycc (allocate) 
call g_diag_ccvv_y60(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_diag_ccvv_y60



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y60(s_c1, i_c1, V2_, Y30_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y30_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y
! Y30(y,y)  <-- 
! (    1.00000000)  V2(c1,y,c1,y) 
do s_y = 0, nir-1
if( &
IEOR(s_y,s_y) == 0 .and. & 
IEOR(s_c1,s_y) == IEOR(s_c1,s_y)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y30_(s_y, s_y)%array(i_y, i_y) =  &
    Y30_(s_y, s_y)%array(i_y, i_y) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_c1, s_y)%array(i_y, i_c1, i_y)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y60



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x60(sa, ia, Y30, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y30(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Ycc(sleft2, Y30, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x60(sa, ia, Ycc, av2_i2, nir, nsym, psym)

deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x60



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x60(s_a, i_a, Y30_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y30_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! Hdiag(w,y,a,a)  <-- 
! (   -2.00000000) Y30(y,y) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_y,s_y) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * Y30_(s_y, s_y)%array(i_y, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x60



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x61(sa, ia, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x61(sa, ia, sw, iw, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x61



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x61(s_a, i_a, s_w, i_w, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (    4.00000000) V2(w,y,w,y) 
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_w,s_y) == IEOR(s_w,s_y)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * V2_(s_y, s_w, s_y)%array(i_y, i_w, i_y)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x61



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x62(sa, ia, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x62(sa, ia, sw, iw, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x62



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x62(s_a, i_a, s_w, i_w, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -2.00000000) V2(w,w,y,y) 
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_w,s_w) == IEOR(s_y,s_y)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * V2_(s_y, s_y, s_w)%array(i_y, i_y, i_w)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x62



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y63(sc1, ic1, V2, Y31, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y31(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y31, nir, nsym, psym) ! -> Ycc (allocate) 
call g_diag_ccvv_y63(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_diag_ccvv_y63



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y63(s_c1, i_c1, V2_, Y31_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y31_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Y31(w,y)  <-- 
! (    1.00000000)  V2(c1,w,c1,y) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == 0 .and. & 
IEOR(s_c1,s_w) == IEOR(s_c1,s_y)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y31_(s_y, s_w)%array(i_y, i_w) =  &
    Y31_(s_y, s_w)%array(i_y, i_w) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_c1, s_w)%array(i_y, i_c1, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y63



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x63(sa, ia, Y31, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y31(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Ycc(sleft2, Y31, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x63(sa, ia, Ycc, av2_i2, nir, nsym, psym)

deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x63



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x63(s_a, i_a, Y31_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y31_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    8.00000000) Y31(w,w) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_w,s_w) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 8.00000000d+00 & 
  * Y31_(s_w, s_w)%array(i_w, i_w)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x63



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y64(sc1, ic1, V2, Y32, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y32(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y32, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ccvv_y64(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ccvv_y64



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y64(s_c1, i_c1, V2_, Y32_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y32_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y32(o1,o2)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y32_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y32_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_o1, s_c1)%array(i_o2, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y64



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x64(sc, ic, Y32, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y32(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y32, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x64(sc, ic, Yaa, av2_i2, d1, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x64



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x64(s_c, i_c, Y32_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y32_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    8.00000000) D1(o1,o2) Y32(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y32_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_ccvv_no0_x64



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y65(sc1, ic1, V2, Y33, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y33(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y33, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ccvv_y65(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ccvv_y65



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y65(s_c1, i_c1, V2_, Y33_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y33_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y33(o1,o2)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y33_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y33_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_o1, s_c1)%array(i_o2, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y65



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x65(sc, ic, Y33, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y33(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y33, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x65(sc, ic, Yaa, av2_i2, d1, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x65



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x65(s_c, i_c, Y33_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y33_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (   -4.00000000) D1(o1,o2) Y33(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y33_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_ccvv_no0_x65



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x66(sc, ic, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x66(sc, ic, sw, iw, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x66



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x66(s_c, i_c, s_w, i_w, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (   -4.00000000) D1(o1,o2) V2(w,w,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_w,s_w) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_w)%array(i_o2, i_o1, i_w)
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

end subroutine g_diag_ccvv_no0_x66



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x67(sc, ic, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x67(sc, ic, sw, iw, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x67



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x67(s_c, i_c, s_w, i_w, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (    4.00000000) D1(o1,o2) V2(w,w,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_w,s_w) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_w)%array(i_o2, i_o1, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x67



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x68(sc, ic, sy, iy, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sy, iy
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sy, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x68(sc, ic, sy, iy, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x68



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x68(s_c, i_c, s_y, i_y, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_y, s_y
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (   -4.00000000) D1(o1,o2) V2(y,y,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_y,s_y) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_y)%array(i_o2, i_o1, i_y)
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

end subroutine g_diag_ccvv_no0_x68



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y69(sc1, ic1, V2, Y34, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y34(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y34, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ccvv_y69(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ccvv_y69



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y69(s_c1, i_c1, V2_, Y34_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y34_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y34(o1,o2)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y34_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y34_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_o1, s_c1)%array(i_o2, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y69



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x69(sa, ia, Y34, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y34(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y34, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x69(sa, ia, Yaa, av2_i2, d1, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x69



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x69(s_a, i_a, Y34_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y34_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    4.00000000) D1(o1,o2) Y34(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y34_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x69



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y70(sc1, ic1, V2, Y35, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y35(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y35, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ccvv_y70(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ccvv_y70



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y70(s_c1, i_c1, V2_, Y35_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y35_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y35(o1,o2)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y35_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y35_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_o1, s_c1)%array(i_o2, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y70



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x70(sa, ia, Y35, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y35(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y35, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x70(sa, ia, Yaa, av2_i2, d1, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x70



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x70(s_a, i_a, Y35_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y35_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -2.00000000) D1(o1,o2) Y35(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y35_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_ccvv_no0_x70



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x71(sa, ia, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x71(sa, ia, sw, iw, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x71



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x71(s_a, i_a, s_w, i_w, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (    1.00000000) D1(o1,o2) V2(w,w,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_w,s_w) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_w)%array(i_o2, i_o1, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x71



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x72(sa, ia, sy, iy, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sy, iy
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sy, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x72(sa, ia, sy, iy, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x72



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x72(s_a, i_a, s_y, i_y, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_y, s_y
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,y,a,a)  <-- 
! (    1.00000000) D1(o1,o2) V2(y,y,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_y,s_y) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_y)%array(i_o2, i_o1, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x72



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x73(sa, ia, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x73(sa, ia, sw, iw, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x73



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x73(s_a, i_a, s_w, i_w, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Hdiag(w,w,a,a)  <-- 
! (   -4.00000000) D1(o1,o2) V2(w,w,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_w,s_w) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_w)%array(i_o2, i_o1, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x73



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x74(sc, ic, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x74(sc, ic, sw, iw, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x74



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x74(s_c, i_c, s_w, i_w, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (   -2.00000000) D1(o1,o2) V2(w,o2,w,o1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_w,s_o2) == IEOR(s_w,s_o1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o1, s_w, s_o2)%array(i_o1, i_w, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x74



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x75(sc, ic, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x75(sc, ic, sw, iw, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x75



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x75(s_c, i_c, s_w, i_w, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    2.00000000) D1(o1,o2) V2(w,o1,w,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_w,s_o1) == IEOR(s_w,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_w, s_o1)%array(i_o2, i_w, i_o1)
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

end subroutine g_diag_ccvv_no0_x75



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x76(sc, ic, sy, iy, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sy, iy
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sy, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x76(sc, ic, sy, iy, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x76



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x76(s_c, i_c, s_y, i_y, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_y, s_y
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    2.00000000) D1(o1,o2) V2(y,o1,y,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_y,s_o1) == IEOR(s_y,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_y, s_o1)%array(i_o2, i_y, i_o1)
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

end subroutine g_diag_ccvv_no0_x76



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y77(sc1, ic1, V2, Y36, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y36(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y36, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ccvv_y77(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ccvv_y77



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y77(s_c1, i_c1, V2_, Y36_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y36_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y36(o1,o2)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y36_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y36_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_c1, s_o1)%array(i_o2, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y77



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x77(sc, ic, Y36, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y36(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y36, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x77(sc, ic, Yaa, av2_i2, d1, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x77



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x77(s_c, i_c, Y36_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y36_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (   -4.00000000) D1(o1,o2) Y36(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y36_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_ccvv_no0_x77



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x78(sa, ia, sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x78(sa, ia, sc, ic, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x78



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x78(s_a, i_a, s_c, i_c, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,c)  <-- 
! (    4.00000000) D1(o1,o2) V2(a,a,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_a) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
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

end subroutine g_diag_ccvv_no0_x78



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x79(sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x79(sc, ic, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x79



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x79(s_c, i_c, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    4.00000000) D1(o1,o2) V2(c,c,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c,s_c) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
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

end subroutine g_diag_ccvv_no0_x79



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x80(sa, ia, sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x80(sa, ia, sc, ic, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x80



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x80(s_a, i_a, s_c, i_c, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,c)  <-- 
! (   -2.00000000) D1(o1,o2) V2(a,o1,o2,a) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_o1) == IEOR(s_o2,s_a)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
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

end subroutine g_diag_ccvv_no0_x80



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x81(sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x81(sc, ic, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x81



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x81(s_c, i_c, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y, s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (   -2.00000000) D1(o1,o2) V2(c,o1,o2,c) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c,s_o1) == IEOR(s_o2,s_c)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
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

end subroutine g_diag_ccvv_no0_x81



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y82(sc1, ic1, V2, Y37, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y37(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y37, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ccvv_y82(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ccvv_y82



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y82(s_c1, i_c1, V2_, Y37_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y37_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y37(o1,o2)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y37_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y37_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_c1, s_o1)%array(i_o2, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y82



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x82(sc, ic, Y37, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y37(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y37, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x82(sc, ic, Yaa, av2_i2, d1, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x82



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x82(s_c, i_c, Y37_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y37_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (    2.00000000) D1(o1,o2) Y37(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y37_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_ccvv_no0_x82



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x83(sa, ia, sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x83(sa, ia, sc, ic, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x83



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x83(s_a, i_a, s_c, i_c, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,w,a,c)  <-- 
! (   -2.00000000) D1(o1,o2) V2(a,a,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_a) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_a)%array(i_o2, i_o1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x83



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x84(sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x84(sc, ic, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x84



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x84(s_c, i_c, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (   -2.00000000) D1(o1,o2) V2(c,c,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c,s_c) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_c)%array(i_o2, i_o1, i_c)
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

end subroutine g_diag_ccvv_no0_x84



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x85(sa, ia, sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x85(sa, ia, sc, ic, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x85



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x85(s_a, i_a, s_c, i_c, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,w,a,c)  <-- 
! (    1.00000000) D1(o1,o2) V2(a,o1,o2,a) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_o1) == IEOR(s_o2,s_a)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x85



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x86(sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x86(sc, ic, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x86



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x86(s_c, i_c, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (    1.00000000) D1(o1,o2) V2(c,o1,o2,c) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c,s_o1) == IEOR(s_o2,s_c)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_c, s_o2, s_o1)%array(i_c, i_o2, i_o1)
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

end subroutine g_diag_ccvv_no0_x86



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x87(sa, ia, sy, iy, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sy, iy
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sy, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x87(sa, ia, sy, iy, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x87



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x87(s_a, i_a, s_y, i_y, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_y, s_y
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,y,a,a)  <-- 
! (   -0.50000000) D1(o1,o2) V2(y,o1,y,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_y,s_o1) == IEOR(s_y,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 0.50000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_y, s_o1)%array(i_o2, i_y, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x87



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x88(sa, ia, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x88(sa, ia, sw, iw, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x88



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x88(s_a, i_a, s_w, i_w, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Hdiag(w,w,a,a)  <-- 
! (    2.00000000) D1(o1,o2) V2(w,o1,w,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_w,s_o1) == IEOR(s_w,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_w, s_o1)%array(i_o2, i_w, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x88



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y89(sc1, ic1, V2, Y38, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y38(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y38, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ccvv_y89(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ccvv_y89



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y89(s_c1, i_c1, V2_, Y38_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y38_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y38(o1,o2)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y38_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y38_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_c1, s_o1)%array(i_o2, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y89



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x89(sa, ia, Y38, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y38(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y38, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x89(sa, ia, Yaa, av2_i2, d1, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x89



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x89(s_a, i_a, Y38_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y38_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (   -2.00000000) D1(o1,o2) Y38(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y38_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x89



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x90(sa, ia, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x90(sa, ia, sw, iw, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x90



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x90(s_a, i_a, s_w, i_w, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -0.50000000) D1(o1,o2) V2(w,o1,w,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_w,s_o1) == IEOR(s_w,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 0.50000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_w, s_o1)%array(i_o2, i_w, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x90



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y91(sc1, ic1, V2, Y39, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y39(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y39, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ccvv_y91(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ccvv_y91



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y91(s_c1, i_c1, V2_, Y39_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y39_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y39(o1,o2)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y39_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y39_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_c1, s_o1)%array(i_o2, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y91



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x91(sa, ia, Y39, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y39(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y39, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x91(sa, ia, Yaa, av2_i2, d1, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x91



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x91(s_a, i_a, Y39_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y39_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (    1.00000000) D1(o1,o2) Y39(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y39_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_ccvv_no0_x91



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y92(sc1, ic1, V2, Y40, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y40(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y40, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ccvv_y92(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y92



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y92(s_c1, i_c1, V2_, Y40_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y40_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_c, i_c
! Y40(a,c)  <-- 
! (    1.00000000)  V2(c1,c1,a,c) 
do s_a = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a,s_c) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_c)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y40_(s_c, s_a)%array(i_c, i_a) =  &
    Y40_(s_c, s_a)%array(i_c, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_a, s_c1)%array(i_c, i_a, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y92



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x92(sa, ia, Y40, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y40(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y40, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x92(sa, ia, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x92



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x92(s_a, i_a, Y40_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y40_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    8.00000000) Y40(a,a) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 8.00000000d+00 & 
  * Y40_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x92



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y93(sc1, ic1, V2, Y41, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y41(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y41, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ccvv_y93(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y93



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y93(s_c1, i_c1, V2_, Y41_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y41_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_c, i_c
! Y41(a,c)  <-- 
! (    1.00000000)  V2(c1,c1,a,c) 
do s_a = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a,s_c) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_c)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y41_(s_c, s_a)%array(i_c, i_a) =  &
    Y41_(s_c, s_a)%array(i_c, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_a, s_c1)%array(i_c, i_a, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y93



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x93(sa, ia, Y41, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y41(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y41, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x93(sa, ia, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x93



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x93(s_a, i_a, Y41_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y41_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -4.00000000) Y41(a,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * Y41_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x93



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x94(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x94(sa, ia, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x94



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x94(s_a, i_a, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (   -8.00000000) V2(a,a,w,w) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == IEOR(s_w,s_w)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 8.00000000d+00 & 
  * V2_(s_w, s_w, s_a)%array(i_w, i_w, i_a)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x94



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x95(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x95(sa, ia, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x95



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x95(s_a, i_a, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (    2.00000000) V2(a,a,w,w) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == IEOR(s_w,s_w)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * V2_(s_w, s_w, s_a)%array(i_w, i_w, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x95



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x96(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x96(sa, ia, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x96



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x96(s_a, i_a, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! Hdiag(w,y,a,a)  <-- 
! (    2.00000000) V2(a,a,y,y) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == IEOR(s_y,s_y)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * V2_(s_y, s_y, s_a)%array(i_y, i_y, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x96



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x97(sa, ia, sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x97(sa, ia, sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x97



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x97(s_a, i_a, s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,c)  <-- 
! (    4.00000000) V2(a,a,w,w) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_a,s_a) == IEOR(s_w,s_w)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 4.00000000d+00 & 
  * V2_(s_w, s_w, s_a)%array(i_w, i_w, i_a)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x97



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x98(sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x98(sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x98



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x98(s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (    4.00000000) V2(c,c,w,w) 
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_c) == IEOR(s_w,s_w)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 4.00000000d+00 & 
  * V2_(s_w, s_w, s_c)%array(i_w, i_w, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x98



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x99(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x99(sa, ia, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x99



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x99(s_a, i_a, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    8.00000000) V2(a,w,w,a) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_w,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 8.00000000d+00 & 
  * V2_(s_a, s_w, s_w)%array(i_a, i_w, i_w)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x99



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x100(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x100(sa, ia, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x100



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x100(s_a, i_a, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -4.00000000) V2(a,w,w,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_w,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * V2_(s_a, s_w, s_w)%array(i_a, i_w, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x100



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x101(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x101(sa, ia, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x101



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x101(s_a, i_a, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! Hdiag(w,y,a,a)  <-- 
! (   -4.00000000) V2(a,y,y,a) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_y) == IEOR(s_y,s_a)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * V2_(s_a, s_y, s_y)%array(i_a, i_y, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x101



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y102(sc1, ic1, V2, Y42, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y42(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y42, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ccvv_y102(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y102



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y102(s_c1, i_c1, V2_, Y42_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y42_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_c, i_c
! Y42(a,c)  <-- 
! (    1.00000000)  V2(c1,a,c1,c) 
do s_a = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a,s_c) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_c)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y42_(s_c, s_a)%array(i_c, i_a) =  &
    Y42_(s_c, s_a)%array(i_c, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_c1, s_a)%array(i_c, i_c1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y102



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x102(sa, ia, Y42, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y42(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y42, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x102(sa, ia, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x102



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x102(s_a, i_a, Y42_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y42_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (    2.00000000) Y42(a,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * Y42_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x102



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x103(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x103(sa, ia, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x103



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x103(s_a, i_a, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    2.00000000) V2(a,w,w,a) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_w,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 2.00000000d+00 & 
  * V2_(s_a, s_w, s_w)%array(i_a, i_w, i_w)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x103



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y104(sc1, ic1, V2, Y43, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y43(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y43, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ccvv_y104(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y104



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y104(s_c1, i_c1, V2_, Y43_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y43_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_c, i_c
! Y43(a,c)  <-- 
! (    1.00000000)  V2(c1,a,c1,c) 
do s_a = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a,s_c) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_c)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y43_(s_c, s_a)%array(i_c, i_a) =  &
    Y43_(s_c, s_a)%array(i_c, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_c1, s_a)%array(i_c, i_c1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y104



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x104(sa, ia, Y43, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y43(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y43, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x104(sa, ia, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x104



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x104(s_a, i_a, Y43_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y43_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (   -4.00000000) Y43(a,a) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 4.00000000d+00 & 
  * Y43_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x104



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x105(sa, ia, sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x105(sa, ia, sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x105



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x105(s_a, i_a, s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,c)  <-- 
! (   -8.00000000) V2(a,w,w,a) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_a,s_w) == IEOR(s_w,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 8.00000000d+00 & 
  * V2_(s_a, s_w, s_w)%array(i_a, i_w, i_w)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x105



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x106(sc, ic, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x106(sc, ic, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x106



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x106(s_c, i_c, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (   -8.00000000) V2(c,w,w,c) 
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_w) == IEOR(s_w,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 8.00000000d+00 & 
  * V2_(s_c, s_w, s_w)%array(i_c, i_w, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x106



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x107(sa, ia, sy, iy, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sy, iy
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sy, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x107(sa, ia, sy, iy, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x107



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x107(s_a, i_a, s_y, i_y, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_y, s_y
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,y,a,a)  <-- 
! (   -0.50000000) D1(o1,o2) V2(y,o1,y,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_y,s_o1) == IEOR(s_y,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 0.50000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_y, s_o1)%array(i_o2, i_y, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x107



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x108(sa, ia, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x108(sa, ia, sw, iw, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x108



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x108(s_a, i_a, s_w, i_w, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Hdiag(w,w,a,a)  <-- 
! (    2.00000000) D1(o1,o2) V2(w,o1,w,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_w,s_o1) == IEOR(s_w,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_w, s_o1)%array(i_o2, i_w, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x108



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y109(sc1, ic1, V2, Y44, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y44(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y44, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ccvv_y109(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ccvv_y109



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y109(s_c1, i_c1, V2_, Y44_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y44_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y44(o1,o2)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y44_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y44_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_c1, s_o1)%array(i_o2, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y109



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x109(sa, ia, Y44, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y44(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y44, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x109(sa, ia, Yaa, av2_i2, d1, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x109



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x109(s_a, i_a, Y44_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y44_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (   -2.00000000) D1(o1,o2) Y44(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y44_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x109



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x110(sa, ia, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x110(sa, ia, sw, iw, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x110



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x110(s_a, i_a, s_w, i_w, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -0.50000000) D1(o1,o2) V2(w,o1,w,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_w,s_o1) == IEOR(s_w,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 0.50000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_w, s_o1)%array(i_o2, i_w, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x110



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y111(sc1, ic1, V2, Y45, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y45(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y45, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ccvv_y111(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ccvv_y111



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y111(s_c1, i_c1, V2_, Y45_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y45_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y45(o1,o2)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y45_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y45_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_c1, s_o1)%array(i_o2, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y111



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x111(sa, ia, Y45, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y45(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y45, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x111(sa, ia, Yaa, av2_i2, d1, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x111



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x111(s_a, i_a, Y45_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y45_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (    1.00000000) D1(o1,o2) Y45(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y45_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_ccvv_no0_x111



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y112(sc1, ic1, V2, Y46, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y46(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y46, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ccvv_y112(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ccvv_y112



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y112(s_c1, i_c1, V2_, Y46_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y46_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y46(o1,o2)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y46_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y46_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_o1, s_c1)%array(i_o2, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y112



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x112(sa, ia, Y46, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y46(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y46, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x112(sa, ia, Yaa, av2_i2, d1, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x112



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x112(s_a, i_a, Y46_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y46_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    4.00000000) D1(o1,o2) Y46(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y46_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x112



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y113(sc1, ic1, V2, Y47, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y47(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y47, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ccvv_y113(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ccvv_y113



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y113(s_c1, i_c1, V2_, Y47_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y47_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Y47(o1,o2)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y47_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y47_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_o1, s_c1)%array(i_o2, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y113



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x113(sa, ia, Y47, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y47(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y47, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x113(sa, ia, Yaa, av2_i2, d1, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x113



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x113(s_a, i_a, Y47_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y47_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -2.00000000) D1(o1,o2) Y47(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y47_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_diag_ccvv_no0_x113



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x114(sa, ia, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x114(sa, ia, sw, iw, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x114



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x114(s_a, i_a, s_w, i_w, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (    1.00000000) D1(o1,o2) V2(w,w,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_w,s_w) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_w)%array(i_o2, i_o1, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x114



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x115(sa, ia, sy, iy, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sy, iy
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sy, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x115(sa, ia, sy, iy, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x115



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x115(s_a, i_a, s_y, i_y, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_y, s_y
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,y,a,a)  <-- 
! (    1.00000000) D1(o1,o2) V2(y,y,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_y,s_y) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_y)%array(i_o2, i_o1, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x115



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x116(sa, ia, sw, iw, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x116(sa, ia, sw, iw, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x116



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x116(s_a, i_a, s_w, i_w, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! Hdiag(w,w,a,a)  <-- 
! (   -4.00000000) D1(o1,o2) V2(w,w,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_w,s_w) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_w)%array(i_o2, i_o1, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x116



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x117(sc, ic, so1, io1, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, so1, io1
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x117(sc, ic, so1, io1, h2_i, av2_i2, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x117



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x117(s_c, i_c, s_o1, i_o1, V2_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_o4, i_o4, s_w, i_w, s_y, i_y
integer :: s_a, i_a
! Hdiag(w,y,a,c)  <-- 
! (    2.00000000) D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * V2_(s_o4, s_o2, s_o3)%array(i_o4, i_o2, i_o3)
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

end subroutine g_diag_ccvv_no0_x117



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x118(sc, ic, so1, io1, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, so1, io1
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x118(sc, ic, so1, io1, h2_i, av2_i2, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x118



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x118(s_c, i_c, s_o1, i_o1, V2_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_o4, i_o4, s_w, i_w, s_a, i_a
! Hdiag(w,w,a,c)  <-- 
! (   -1.00000000) D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
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

end subroutine g_diag_ccvv_no0_x118



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x119(sa, ia, so1, io1, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, so1, io1
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x119(sa, ia, so1, io1, h2_i, av2_i2, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x119



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x119(s_a, i_a, s_o1, i_o1, V2_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_o4, i_o4, s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    2.00000000) D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 2.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
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

end subroutine g_diag_ccvv_no0_x119



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x120(sa, ia, so1, io1, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, so1, io1
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x120(sa, ia, so1, io1, h2_i, av2_i2, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x120



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x120(s_a, i_a, s_o1, i_o1, V2_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_o4, i_o4, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -1.00000000) D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
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

end subroutine g_diag_ccvv_no0_x120



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x121(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x121(sa, ia, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x121



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x121(s_a, i_a, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    4.00000000) D1(o1,o2) V2(a,a,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_a) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_a)%array(i_o2, i_o1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x121



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x122(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x122(sa, ia, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x122



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x122(s_a, i_a, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -2.00000000) D1(o1,o2) V2(a,a,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_a) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
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

end subroutine g_diag_ccvv_no0_x122



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x123(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x123(sa, ia, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x123



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x123(s_a, i_a, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (   -2.00000000) D1(o1,o2) V2(a,o1,o2,a) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_o1) == IEOR(s_o2,s_a)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x123



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x124(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x124(sa, ia, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x124



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x124(s_a, i_a, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (    1.00000000) D1(o1,o2) V2(a,o1,o2,a) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_o1) == IEOR(s_o2,s_a)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
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

end subroutine g_diag_ccvv_no0_x124



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x125(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x125(sa, ia, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x125



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x125(s_a, i_a, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    8.00000000) V2(a,w,w,a) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_w,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 8.00000000d+00 & 
  * V2_(s_a, s_w, s_w)%array(i_a, i_w, i_w)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x125



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x126(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x126(sa, ia, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x126



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x126(s_a, i_a, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -4.00000000) V2(a,w,w,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_w,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * V2_(s_a, s_w, s_w)%array(i_a, i_w, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x126



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x127(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x127(sa, ia, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x127



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x127(s_a, i_a, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! Hdiag(w,y,a,a)  <-- 
! (   -4.00000000) V2(a,y,y,a) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_y) == IEOR(s_y,s_a)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * V2_(s_a, s_y, s_y)%array(i_a, i_y, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x127



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y128(sc1, ic1, V2, Y48, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y48(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y48, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ccvv_y128(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y128



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y128(s_c1, i_c1, V2_, Y48_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y48_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_c, i_c
! Y48(a,c)  <-- 
! (    1.00000000)  V2(c1,a,c1,c) 
do s_a = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a,s_c) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_c)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y48_(s_c, s_a)%array(i_c, i_a) =  &
    Y48_(s_c, s_a)%array(i_c, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_c1, s_a)%array(i_c, i_c1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y128



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x128(sa, ia, Y48, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y48(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y48, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x128(sa, ia, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x128



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x128(s_a, i_a, Y48_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y48_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (    2.00000000) Y48(a,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * Y48_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x128



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x129(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x129(sa, ia, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x129



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x129(s_a, i_a, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    2.00000000) V2(a,w,w,a) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_w,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 2.00000000d+00 & 
  * V2_(s_a, s_w, s_w)%array(i_a, i_w, i_w)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x129



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y130(sc1, ic1, V2, Y49, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y49(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y49, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ccvv_y130(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y130



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y130(s_c1, i_c1, V2_, Y49_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y49_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_c, i_c
! Y49(a,c)  <-- 
! (    1.00000000)  V2(c1,a,c1,c) 
do s_a = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a,s_c) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_c)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y49_(s_c, s_a)%array(i_c, i_a) =  &
    Y49_(s_c, s_a)%array(i_c, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_c1, s_a)%array(i_c, i_c1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y130



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x130(sa, ia, Y49, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y49(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y49, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x130(sa, ia, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x130



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x130(s_a, i_a, Y49_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y49_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (   -4.00000000) Y49(a,a) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 4.00000000d+00 & 
  * Y49_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x130



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y131(sc1, ic1, V2, Y50, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y50(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y50, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ccvv_y131(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y131



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y131(s_c1, i_c1, V2_, Y50_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y50_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_c, i_c
! Y50(a,c)  <-- 
! (    1.00000000)  V2(c1,c1,a,c) 
do s_a = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a,s_c) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_c)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y50_(s_c, s_a)%array(i_c, i_a) =  &
    Y50_(s_c, s_a)%array(i_c, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_a, s_c1)%array(i_c, i_a, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y131



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x131(sa, ia, Y50, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y50(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y50, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x131(sa, ia, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x131



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x131(s_a, i_a, Y50_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y50_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    8.00000000) Y50(a,a) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 8.00000000d+00 & 
  * Y50_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x131



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_y132(sc1, ic1, V2, Y51, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y51(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y51, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ccvv_y132(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ccvv_y132



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_y132(s_c1, i_c1, V2_, Y51_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y51_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_c, i_c
! Y51(a,c)  <-- 
! (    1.00000000)  V2(c1,c1,a,c) 
do s_a = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a,s_c) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_c)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
Y51_(s_c, s_a)%array(i_c, i_a) =  &
    Y51_(s_c, s_a)%array(i_c, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_c, s_a, s_c1)%array(i_c, i_a, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_y132



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x132(sa, ia, Y51, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y51(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y51, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ccvv_no0_x132(sa, ia, Yvv, av2_i2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x132



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x132(s_a, i_a, Y51_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y51_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -4.00000000) Y51(a,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * Y51_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x132



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x133(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x133(sa, ia, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x133



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x133(s_a, i_a, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (   -8.00000000) V2(a,a,w,w) 
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == IEOR(s_w,s_w)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 8.00000000d+00 & 
  * V2_(s_w, s_w, s_a)%array(i_w, i_w, i_a)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x133



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x134(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x134(sa, ia, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x134



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x134(s_a, i_a, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (    2.00000000) V2(a,a,w,w) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == IEOR(s_w,s_w)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * V2_(s_w, s_w, s_a)%array(i_w, i_w, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x134



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x135(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x135(sa, ia, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x135



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x135(s_a, i_a, V2_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! Hdiag(w,y,a,a)  <-- 
! (    2.00000000) V2(a,a,y,y) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_a,s_a) == IEOR(s_y,s_y)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * V2_(s_y, s_y, s_a)%array(i_y, i_y, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x135



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x136(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x136(sa, ia, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x136



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x136(s_a, i_a, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (   -2.00000000) D1(o1,o2) V2(a,o1,o2,a) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_o1) == IEOR(s_o2,s_a)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x136



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x137(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x137(sa, ia, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x137



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x137(s_a, i_a, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (    1.00000000) D1(o1,o2) V2(a,o1,o2,a) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_o1) == IEOR(s_o2,s_a)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
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

end subroutine g_diag_ccvv_no0_x137



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x138(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x138(sa, ia, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x138



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x138(s_a, i_a, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! Hdiag(w,w,a,a)  <-- 
! (    4.00000000) D1(o1,o2) V2(a,a,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_w) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_a) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) =  &
    Hdiag_(s_a, s_w, s_w)%array(i_a, i_w, i_w) &
  + 4.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_a)%array(i_o2, i_o1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ccvv_no0_x138



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ccvv_no0_x139(sa, ia, V2, Hdiag, nir, nsym, psym)

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
call g_diag_ccvv_no0_x139(sa, ia, h2_i, av2_i2, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ccvv_no0_x139



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ccvv_no0_x139(s_a, i_a, V2_, Hdiag_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w, s_y, i_y
! Hdiag(w,y,a,a)  <-- 
! (   -2.00000000) D1(o1,o2) V2(a,a,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_a) .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_a) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    Hdiag_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
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

end subroutine g_diag_ccvv_no0_x139

