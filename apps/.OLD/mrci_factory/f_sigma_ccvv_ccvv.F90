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
subroutine g_if_sigma_ccvv_ccvv_y0(h, Y0, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y0(h1, Y0, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_sigma_ccvv_ccvv_y0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y0(h_, Y0_, nir, nsym, psym)

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

end subroutine g_sigma_ccvv_ccvv_y0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x0(sc, ic, Y0, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y0
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x0(sc, ic, Y0, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x0(s_c, i_c, Y0, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y0
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    8.00000000) Y0 T2(w,y,a,c) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_a,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * Y0 & 
  * T2_(s_a, s_y, s_w)%array(i_a, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y1(h, Y1, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y1(h1, Y1, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_sigma_ccvv_ccvv_y1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y1(h_, Y1_, nir, nsym, psym)

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

end subroutine g_sigma_ccvv_ccvv_y1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x1(sc, ic, Y1, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y1
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x1(sc, ic, Y1, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x1(s_c, i_c, Y1, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y1
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) Y1 T2(y,w,a,c) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_a,s_c)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * Y1 & 
  * T2_(s_a, s_w, s_y)%array(i_a, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x2(sc, ic, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x2(sc, ic, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x2(s_c, i_c, T2_, h_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_y, i_y, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(c1,y,a,c) h(c1,w) 
do s_c1 = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_y) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_w) == 0) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_a, s_y, s_c1)%array(i_a, i_y, i_c1) & 
  * h_(s_w, s_c1)%array(i_w, i_c1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x3(sc, ic, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x3(sc, ic, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x3(s_c, i_c, T2_, h_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_c1, i_c1, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(y,c1,a,c) h(c1,w) 
do s_y = 0, nir-1
do s_c1 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c1) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_w) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_a, s_c1, s_y)%array(i_a, i_c1, i_y) & 
  * h_(s_w, s_c1)%array(i_w, i_c1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x4(sc, ic, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x4(sc, ic, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x4(s_c, i_c, T2_, h_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_w, i_w, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(c1,w,a,c) h(c1,y) 
do s_c1 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_w) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_y) == 0) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_a, s_w, s_c1)%array(i_a, i_w, i_c1) & 
  * h_(s_y, s_c1)%array(i_y, i_c1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x5(sc, ic, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x5(sc, ic, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x5(s_c, i_c, T2_, h_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c1, i_c1, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(w,c1,a,c) h(c1,y) 
do s_w = 0, nir-1
do s_c1 = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c1) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_y) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_a, s_c1, s_w)%array(i_a, i_c1, i_w) & 
  * h_(s_y, s_c1)%array(i_y, i_c1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y6(h, Y2, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y6(h1, Y2, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_sigma_ccvv_ccvv_y6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y6(h_, Y2_, nir, nsym, psym)

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

end subroutine g_sigma_ccvv_ccvv_y6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x6(sa, ia, sc, ic, Y2, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: Y2
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x6(sa, ia, sc, ic, Y2, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x6(s_a, i_a, s_c, i_c, Y2, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y2
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    8.00000000) Y2 T2(y,w,c,a) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_c,s_a)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * Y2 & 
  * T2_(s_c, s_w, s_y)%array(i_c, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y7(h, Y3, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y7(h1, Y3, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_sigma_ccvv_ccvv_y7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y7(h_, Y3_, nir, nsym, psym)

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

end subroutine g_sigma_ccvv_ccvv_y7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x7(sa, ia, sc, ic, Y3, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: Y3
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x7(sa, ia, sc, ic, Y3, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x7(s_a, i_a, s_c, i_c, Y3, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y3
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -4.00000000) Y3 T2(w,y,c,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_c,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * Y3 & 
  * T2_(s_c, s_y, s_w)%array(i_c, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x8(sa, ia, sc, ic, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x8(sa, ia, sc, ic, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x8(s_a, i_a, s_c, i_c, T2_, h_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_c1, i_c1, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(y,c1,c,a) h(c1,w) 
do s_y = 0, nir-1
do s_c1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c1) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_w) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_c, s_c1, s_y)%array(i_c, i_c1, i_y) & 
  * h_(s_w, s_c1)%array(i_w, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x9(sa, ia, sc, ic, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x9(sa, ia, sc, ic, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x9(s_a, i_a, s_c, i_c, T2_, h_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(c1,y,c,a) h(c1,w) 
do s_c1 = 0, nir-1
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_y) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_w) == 0) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_c, s_y, s_c1)%array(i_c, i_y, i_c1) & 
  * h_(s_w, s_c1)%array(i_w, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x10(sa, ia, sc, ic, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x10(sa, ia, sc, ic, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x10(s_a, i_a, s_c, i_c, T2_, h_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c1, i_c1, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(w,c1,c,a) h(c1,y) 
do s_w = 0, nir-1
do s_c1 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c1) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_y) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_c, s_c1, s_w)%array(i_c, i_c1, i_w) & 
  * h_(s_y, s_c1)%array(i_y, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x11(sa, ia, sc, ic, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x11(sa, ia, sc, ic, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x11(s_a, i_a, s_c, i_c, T2_, h_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(c1,w,c,a) h(c1,y) 
do s_c1 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_w) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_y) == 0) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_c, s_w, s_c1)%array(i_c, i_w, i_c1) & 
  * h_(s_y, s_c1)%array(i_y, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x12(h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: h(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = 0
call g_sigma_ccvv_ccvv_no0_x12(h1, x, d1, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_sigma_ccvv_ccvv_no0_x12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x12(h_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X()  <-- 
! (    1.00000000)  D1(o1,o2) h(o2,o1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o2,s_o1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * h_(s_o1, s_o2)%array(i_o1, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x12(sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x12(sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x12(s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    4.00000000) X T2(w,y,a,c) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_a,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * X & 
  * T2_(s_a, s_y, s_w)%array(i_a, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x13(h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: h(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = 0
call g_sigma_ccvv_ccvv_no0_x13(h1, x, d1, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_sigma_ccvv_ccvv_no0_x13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x13(h_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X()  <-- 
! (    1.00000000)  D1(o1,o2) h(o2,o1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o2,s_o1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * h_(s_o1, s_o2)%array(i_o1, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x13(sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x13(sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x13(s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -2.00000000) X T2(y,w,a,c) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_a,s_c)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * X & 
  * T2_(s_a, s_w, s_y)%array(i_a, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y14(sc1, ic1, V2, Y4, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y14(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_y14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y14(s_c1, i_c1, V2_, Y4_, nir, nsym, psym)

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

end subroutine g_sigma_ccvv_ccvv_y14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x14(Y4, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: Y4(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y4, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = 0
call g_sigma_ccvv_ccvv_no0_x14(Yaa, x, d1, nir, nsym, psym)

deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_no0_x14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x14(Y4_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y4_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X()  <-- 
! (    1.00000000)  D1(o1,o2) Y4(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y4_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x14(sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x14(sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x14(s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) X T2(y,w,a,c) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_a,s_c)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * X & 
  * T2_(s_a, s_w, s_y)%array(i_a, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y15(sc1, ic1, V2, Y5, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y15(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_y15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y15(s_c1, i_c1, V2_, Y5_, nir, nsym, psym)

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

end subroutine g_sigma_ccvv_ccvv_y15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x15(Y5, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: Y5(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y5, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = 0
call g_sigma_ccvv_ccvv_no0_x15(Yaa, x, d1, nir, nsym, psym)

deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_no0_x15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x15(Y5_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y5_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X()  <-- 
! (    1.00000000)  D1(o1,o2) Y5(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y5_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x15(sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x15(sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x15(s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    2.00000000) X T2(y,w,a,c) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_a,s_c)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * X & 
  * T2_(s_a, s_w, s_y)%array(i_a, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x16(h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: h(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = 0
call g_sigma_ccvv_ccvv_no0_x16(h1, x, d1, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_sigma_ccvv_ccvv_no0_x16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x16(h_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X()  <-- 
! (    1.00000000)  D1(o1,o2) h(o2,o1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o2,s_o1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * h_(s_o1, s_o2)%array(i_o1, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x16(sa, ia, sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x16(sa, ia, sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x16(s_a, i_a, s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    4.00000000) X T2(y,w,c,a) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_c,s_a)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * X & 
  * T2_(s_c, s_w, s_y)%array(i_c, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x17(h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: h(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = 0
call g_sigma_ccvv_ccvv_no0_x17(h1, x, d1, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_sigma_ccvv_ccvv_no0_x17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x17(h_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X()  <-- 
! (    1.00000000)  D1(o1,o2) h(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * h_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x17(sa, ia, sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x17(sa, ia, sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x17(s_a, i_a, s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -2.00000000) X T2(w,y,c,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_c,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * X & 
  * T2_(s_c, s_y, s_w)%array(i_c, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y18(sc1, ic1, V2, Y6, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y18(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_y18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y18(s_c1, i_c1, V2_, Y6_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o2, i_o2
! Y6(o1,o2)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y6_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y6_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_o1, s_c1)%array(i_o2, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x18(Y6, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: Y6(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y6, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = 0
call g_sigma_ccvv_ccvv_no0_x18(Yaa, x, d1, nir, nsym, psym)

deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_no0_x18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x18(Y6_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y6_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X()  <-- 
! (    1.00000000)  D1(o1,o2) Y6(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y6_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x18(sa, ia, sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x18(sa, ia, sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x18(s_a, i_a, s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -4.00000000) X T2(w,y,c,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_c,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * X & 
  * T2_(s_c, s_y, s_w)%array(i_c, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y19(sc1, ic1, V2, Y7, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y19(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_y19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y19(s_c1, i_c1, V2_, Y7_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o2, i_o2
! Y7(o1,o2)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y7_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y7_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_c1, s_o1)%array(i_o2, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x19(Y7, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: Y7(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y7, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = 0
call g_sigma_ccvv_ccvv_no0_x19(Yaa, x, d1, nir, nsym, psym)

deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_no0_x19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x19(Y7_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y7_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X()  <-- 
! (    1.00000000)  D1(o1,o2) Y7(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y7_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x19(sa, ia, sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x19(sa, ia, sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x19(s_a, i_a, s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    2.00000000) X T2(w,y,c,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_c,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * X & 
  * T2_(s_c, s_y, s_w)%array(i_c, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x20(sa, ia, sc, ic, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x20(sa, ia, sc, ic, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x20(s_a, i_a, s_c, i_c, T2_, h_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_v1, i_v1
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(y,w,v1,a) h(c,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_v1, s_w, s_y)%array(i_v1, i_w, i_y) & 
  * h_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x21(sa, ia, sc, ic, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x21(sa, ia, sc, ic, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x21(s_a, i_a, s_c, i_c, T2_, h_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_v1, i_v1
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(w,y,v1,a) h(c,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_v1, s_y, s_w)%array(i_v1, i_y, i_w) & 
  * h_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y22(sc1, ic1, V2, Y8, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y22(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y22



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y22(s_c1, i_c1, V2_, Y8_, nir, nsym, psym)

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

integer :: s_c, i_c, s_v1, i_v1
! Y8(c,v1)  <-- 
! (    1.00000000)  V2(c1,c1,c,v1) 
do s_c = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c,s_v1)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y8_(s_v1, s_c)%array(i_v1, i_c) =  &
    Y8_(s_v1, s_c)%array(i_v1, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c, s_c1)%array(i_v1, i_c, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y22



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x22(sa, ia, sc, ic, T2, Y8, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), Y8(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y8, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x22(sa, ia, sc, ic, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x22



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x22(s_a, i_a, s_c, i_c, T2_, Y8_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: Y8_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_v1, i_v1
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(w,y,v1,a) Y8(c,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_v1, s_y, s_w)%array(i_v1, i_y, i_w) & 
  * Y8_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x22



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y23(sc1, ic1, V2, Y9, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y23(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y23



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y23(s_c1, i_c1, V2_, Y9_, nir, nsym, psym)

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

integer :: s_c, i_c, s_v1, i_v1
! Y9(c,v1)  <-- 
! (    1.00000000)  V2(c1,c,c1,v1) 
do s_c = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_c1,s_c) == IEOR(s_c1,s_v1)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y9_(s_v1, s_c)%array(i_v1, i_c) =  &
    Y9_(s_v1, s_c)%array(i_v1, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_c)%array(i_v1, i_c1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y23



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x23(sa, ia, sc, ic, T2, Y9, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), Y9(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y9, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x23(sa, ia, sc, ic, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x23



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x23(s_a, i_a, s_c, i_c, T2_, Y9_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: Y9_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_v1, i_v1
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(w,y,v1,a) Y9(c,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_v1, s_y, s_w)%array(i_v1, i_y, i_w) & 
  * Y9_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x23



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x24(sc, ic, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x24(sc, ic, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x24



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x24(s_c, i_c, T2_, h_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_v1, i_v1, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(w,y,v1,c) h(a,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_v1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_v1, s_y, s_w)%array(i_v1, i_y, i_w) & 
  * h_(s_v1, s_a)%array(i_v1, i_a)
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

end subroutine g_sigma_ccvv_ccvv_no0_x24



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x25(sc, ic, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x25(sc, ic, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x25



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x25(s_c, i_c, T2_, h_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_v1, i_v1, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(y,w,v1,c) h(a,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_v1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_v1, s_w, s_y)%array(i_v1, i_w, i_y) & 
  * h_(s_v1, s_a)%array(i_v1, i_a)
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

end subroutine g_sigma_ccvv_ccvv_no0_x25



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y26(sc1, ic1, V2, Y10, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y10, nir, nsym, psym) ! -> Yvv (allocate) 
call g_sigma_ccvv_ccvv_y26(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y26



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y26(s_c1, i_c1, V2_, Y10_, nir, nsym, psym)

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

integer :: s_a, i_a, s_v1, i_v1
! Y10(a,v1)  <-- 
! (    1.00000000)  V2(c1,c1,a,v1) 
do s_a = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_v1)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y10_(s_v1, s_a)%array(i_v1, i_a) =  &
    Y10_(s_v1, s_a)%array(i_v1, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_a, s_c1)%array(i_v1, i_a, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y26



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x26(sc, ic, T2, Y10, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), Y10(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y10, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x26(sc, ic, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x26



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x26(s_c, i_c, T2_, Y10_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y10_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_v1, i_v1, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(y,w,v1,c) Y10(a,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_v1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_v1, s_w, s_y)%array(i_v1, i_w, i_y) & 
  * Y10_(s_v1, s_a)%array(i_v1, i_a)
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

end subroutine g_sigma_ccvv_ccvv_no0_x26



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y27(sc1, ic1, V2, Y11, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y11, nir, nsym, psym) ! -> Yvv (allocate) 
call g_sigma_ccvv_ccvv_y27(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y27



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y27(s_c1, i_c1, V2_, Y11_, nir, nsym, psym)

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

integer :: s_a, i_a, s_v1, i_v1
! Y11(a,v1)  <-- 
! (    1.00000000)  V2(c1,a,c1,v1) 
do s_a = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_v1)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y11_(s_v1, s_a)%array(i_v1, i_a) =  &
    Y11_(s_v1, s_a)%array(i_v1, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_a)%array(i_v1, i_c1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y27



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x27(sc, ic, T2, Y11, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), Y11(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y11, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x27(sc, ic, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x27



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x27(s_c, i_c, T2_, Y11_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y11_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_v1, i_v1, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(y,w,v1,c) Y11(a,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_v1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_v1, s_w, s_y)%array(i_v1, i_w, i_y) & 
  * Y11_(s_v1, s_a)%array(i_v1, i_a)
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

end subroutine g_sigma_ccvv_ccvv_no0_x27



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x28(sc, ic, sv1, iv1, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x28(sc, ic, sv1, iv1, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x28



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x28(s_c, i_c, s_v1, i_v1, T2_, h_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(w,y,a,v1) h(c,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) & 
  * h_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x28



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x29(sc, ic, sv1, iv1, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x29(sc, ic, sv1, iv1, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x29



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x29(s_c, i_c, s_v1, i_v1, T2_, h_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(y,w,a,v1) h(c,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_a, s_w, s_y)%array(i_a, i_w, i_y) & 
  * h_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x29



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y30(sc1, ic1, V2, Y12, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y30(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y30



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y30(s_c1, i_c1, V2_, Y12_, nir, nsym, psym)

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

integer :: s_c, i_c, s_v1, i_v1
! Y12(c,v1)  <-- 
! (    1.00000000)  V2(c1,c1,c,v1) 
do s_c = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c,s_v1)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y12_(s_v1, s_c)%array(i_v1, i_c) =  &
    Y12_(s_v1, s_c)%array(i_v1, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c, s_c1)%array(i_v1, i_c, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y30



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x30(sc, ic, sv1, iv1, T2, Y12, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y12(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y12, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x30(sc, ic, sv1, iv1, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x30



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x30(s_c, i_c, s_v1, i_v1, T2_, Y12_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y12_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(y,w,a,v1) Y12(c,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_a, s_w, s_y)%array(i_a, i_w, i_y) & 
  * Y12_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x30



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y31(sc1, ic1, V2, Y13, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y31(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y31



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y31(s_c1, i_c1, V2_, Y13_, nir, nsym, psym)

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

integer :: s_c, i_c, s_v1, i_v1
! Y13(c,v1)  <-- 
! (    1.00000000)  V2(c1,c,c1,v1) 
do s_c = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_c1,s_c) == IEOR(s_c1,s_v1)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y13_(s_v1, s_c)%array(i_v1, i_c) =  &
    Y13_(s_v1, s_c)%array(i_v1, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_c)%array(i_v1, i_c1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y31



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x31(sc, ic, sv1, iv1, T2, Y13, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y13(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y13, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x31(sc, ic, sv1, iv1, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x31



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x31(s_c, i_c, s_v1, i_v1, T2_, Y13_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y13_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(y,w,a,v1) Y13(c,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_a, s_w, s_y)%array(i_a, i_w, i_y) & 
  * Y13_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x31



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x32(sc, ic, sv1, iv1, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x32(sc, ic, sv1, iv1, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x32



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x32(s_c, i_c, s_v1, i_v1, T2_, h_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(y,w,c,v1) h(a,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_c,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_c, s_w, s_y)%array(i_c, i_w, i_y) & 
  * h_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x32



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x33(sc, ic, sv1, iv1, T2, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x33(sc, ic, sv1, iv1, av2_i, h1, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x33



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x33(s_c, i_c, s_v1, i_v1, T2_, h_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(w,y,c,v1) h(a,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_c,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_c, s_y, s_w)%array(i_c, i_y, i_w) & 
  * h_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x33



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y34(sc1, ic1, V2, Y14, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y34(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y34



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y34(s_c1, i_c1, V2_, Y14_, nir, nsym, psym)

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

end subroutine g_sigma_ccvv_ccvv_y34



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x34(sc, ic, sv1, iv1, T2, Y14, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y14(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y14, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x34(sc, ic, sv1, iv1, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x34



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x34(s_c, i_c, s_v1, i_v1, T2_, Y14_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y14_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(w,y,c,v1) Y14(a,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_c,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_c, s_y, s_w)%array(i_c, i_y, i_w) & 
  * Y14_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x34



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y35(sc1, ic1, V2, Y15, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y35(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y35



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y35(s_c1, i_c1, V2_, Y15_, nir, nsym, psym)

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

end subroutine g_sigma_ccvv_ccvv_y35



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x35(sc, ic, sv1, iv1, T2, Y15, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y15(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y15, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x35(sc, ic, sv1, iv1, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x35



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x35(s_c, i_c, s_v1, i_v1, T2_, Y15_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y15_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(w,y,c,v1) Y15(a,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_c,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_c, s_y, s_w)%array(i_c, i_y, i_w) & 
  * Y15_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x35



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y36(sc1, ic1, V2, Y16, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y16
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_sigma_ccvv_ccvv_y36(sc1, ic1, h2_i, Y16, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_y36



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y36(s_c1, i_c1, V2_, Y16_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y16_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y16()  <-- 
! (    1.00000000)  V2(c1,c1,c2,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == IEOR(s_c2,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y16_ = Y16_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c2, s_c1)%array(i_c2, i_c2, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y36



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x36(sc, ic, Y16, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y16
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x36(sc, ic, Y16, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x36



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x36(s_c, i_c, Y16, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y16
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    8.00000000) Y16 T2(w,y,a,c) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_a,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * Y16 & 
  * T2_(s_a, s_y, s_w)%array(i_a, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x36



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y37(sc1, ic1, V2, Y17, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y17
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_sigma_ccvv_ccvv_y37(sc1, ic1, h2_i, Y17, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_y37



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y37(s_c1, i_c1, V2_, Y17_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y17_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y17()  <-- 
! (    1.00000000)  V2(c1,c2,c1,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c2) == IEOR(s_c1,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y17_ = Y17_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c1, s_c2)%array(i_c2, i_c1, i_c2)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y37



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x37(sc, ic, Y17, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y17
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x37(sc, ic, Y17, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x37



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x37(s_c, i_c, Y17, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y17
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) Y17 T2(w,y,a,c) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_a,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * Y17 & 
  * T2_(s_a, s_y, s_w)%array(i_a, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x37



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y38(sc1, ic1, V2, Y18, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y18, nir, nsym, psym) ! -> Ycc (allocate) 
call g_sigma_ccvv_ccvv_y38(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y38



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y38(s_c1, i_c1, V2_, Y18_, nir, nsym, psym)

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

integer :: s_c2, i_c2, s_y, i_y
! Y18(c2,y)  <-- 
! (    1.00000000)  V2(c1,c1,c2,y) 
do s_c2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c2,s_y) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c2,s_y)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y18_(s_y, s_c2)%array(i_y, i_c2) =  &
    Y18_(s_y, s_c2)%array(i_y, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_c2, s_c1)%array(i_y, i_c2, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y38



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x38(sc, ic, T2, Y18, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), Y18(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y18, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x38(sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x38



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x38(s_c, i_c, T2_, Y18_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y18_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c2, i_c2, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -8.00000000) T2(w,c2,a,c) Y18(c2,y) 
do s_w = 0, nir-1
do s_c2 = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c2) == IEOR(s_a,s_c) .and. &
IEOR(s_c2,s_y) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 8.00000000d+00 & 
  * T2_(s_a, s_c2, s_w)%array(i_a, i_c2, i_w) & 
  * Y18_(s_y, s_c2)%array(i_y, i_c2)
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

end subroutine g_sigma_ccvv_ccvv_no0_x38



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y39(sc1, ic1, V2, Y19, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y19, nir, nsym, psym) ! -> Ycc (allocate) 
call g_sigma_ccvv_ccvv_y39(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y39



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y39(s_c1, i_c1, V2_, Y19_, nir, nsym, psym)

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

integer :: s_c2, i_c2, s_y, i_y
! Y19(c2,y)  <-- 
! (    1.00000000)  V2(c1,c2,c1,y) 
do s_c2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c2,s_y) == 0 .and. & 
IEOR(s_c1,s_c2) == IEOR(s_c1,s_y)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y19_(s_y, s_c2)%array(i_y, i_c2) =  &
    Y19_(s_y, s_c2)%array(i_y, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_c1, s_c2)%array(i_y, i_c1, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y39



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x39(sc, ic, T2, Y19, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), Y19(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y19, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x39(sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x39



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x39(s_c, i_c, T2_, Y19_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y19_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c2, i_c2, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(w,c2,a,c) Y19(c2,y) 
do s_w = 0, nir-1
do s_c2 = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c2) == IEOR(s_a,s_c) .and. &
IEOR(s_c2,s_y) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_a, s_c2, s_w)%array(i_a, i_c2, i_w) & 
  * Y19_(s_y, s_c2)%array(i_y, i_c2)
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

end subroutine g_sigma_ccvv_ccvv_no0_x39



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y40(sc1, ic1, V2, Y20, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y20
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_sigma_ccvv_ccvv_y40(sc1, ic1, h2_i, Y20, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_y40



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y40(s_c1, i_c1, V2_, Y20_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y20_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y20()  <-- 
! (    1.00000000)  V2(c1,c1,c2,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == IEOR(s_c2,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y20_ = Y20_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c2, s_c1)%array(i_c2, i_c2, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y40



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x40(sc, ic, Y20, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y20
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x40(sc, ic, Y20, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x40



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x40(s_c, i_c, Y20, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y20
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) Y20 T2(y,w,a,c) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_a,s_c)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * Y20 & 
  * T2_(s_a, s_w, s_y)%array(i_a, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x40



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y41(sc1, ic1, V2, Y21, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y21
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_sigma_ccvv_ccvv_y41(sc1, ic1, h2_i, Y21, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_y41



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y41(s_c1, i_c1, V2_, Y21_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y21_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y21()  <-- 
! (    1.00000000)  V2(c1,c2,c1,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c2) == IEOR(s_c1,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y21_ = Y21_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c1, s_c2)%array(i_c2, i_c1, i_c2)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y41



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x41(sc, ic, Y21, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y21
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x41(sc, ic, Y21, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x41



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x41(s_c, i_c, Y21, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y21
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    2.00000000) Y21 T2(y,w,a,c) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_a,s_c)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * Y21 & 
  * T2_(s_a, s_w, s_y)%array(i_a, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x41



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y42(sc1, ic1, V2, Y22, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y22(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y22, nir, nsym, psym) ! -> Ycc (allocate) 
call g_sigma_ccvv_ccvv_y42(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y42



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y42(s_c1, i_c1, V2_, Y22_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y22_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_y, i_y
! Y22(c2,y)  <-- 
! (    1.00000000)  V2(c1,c1,c2,y) 
do s_c2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c2,s_y) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c2,s_y)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y22_(s_y, s_c2)%array(i_y, i_c2) =  &
    Y22_(s_y, s_c2)%array(i_y, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_c2, s_c1)%array(i_y, i_c2, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y42



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x42(sc, ic, T2, Y22, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), Y22(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y22, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x42(sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x42



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x42(s_c, i_c, T2_, Y22_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y22_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_w, i_w, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(c2,w,a,c) Y22(c2,y) 
do s_c2 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c2,s_w) == IEOR(s_a,s_c) .and. &
IEOR(s_c2,s_y) == 0) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_a, s_w, s_c2)%array(i_a, i_w, i_c2) & 
  * Y22_(s_y, s_c2)%array(i_y, i_c2)
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

end subroutine g_sigma_ccvv_ccvv_no0_x42



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y43(sc1, ic1, V2, Y23, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y23(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y23, nir, nsym, psym) ! -> Ycc (allocate) 
call g_sigma_ccvv_ccvv_y43(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y43



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y43(s_c1, i_c1, V2_, Y23_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y23_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_y, i_y
! Y23(c2,y)  <-- 
! (    1.00000000)  V2(c1,c2,c1,y) 
do s_c2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c2,s_y) == 0 .and. & 
IEOR(s_c1,s_c2) == IEOR(s_c1,s_y)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y23_(s_y, s_c2)%array(i_y, i_c2) =  &
    Y23_(s_y, s_c2)%array(i_y, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_c1, s_c2)%array(i_y, i_c1, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y43



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x43(sc, ic, T2, Y23, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), Y23(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y23, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x43(sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x43



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x43(s_c, i_c, T2_, Y23_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y23_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_w, i_w, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(c2,w,a,c) Y23(c2,y) 
do s_c2 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c2,s_w) == IEOR(s_a,s_c) .and. &
IEOR(s_c2,s_y) == 0) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_a, s_w, s_c2)%array(i_a, i_w, i_c2) & 
  * Y23_(s_y, s_c2)%array(i_y, i_c2)
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

end subroutine g_sigma_ccvv_ccvv_no0_x43



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y44(sc1, ic1, V2, Y24, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y44(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y44



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y44(s_c1, i_c1, V2_, Y24_, nir, nsym, psym)

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

integer :: s_c2, i_c2, s_w, i_w
! Y24(c2,w)  <-- 
! (    1.00000000)  V2(c1,c1,c2,w) 
do s_c2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c2,s_w) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c2,s_w)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Y24_(s_w, s_c2)%array(i_w, i_c2) =  &
    Y24_(s_w, s_c2)%array(i_w, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_w, s_c2, s_c1)%array(i_w, i_c2, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y44



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x44(sc, ic, T2, Y24, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), Y24(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y24, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x44(sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x44



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x44(s_c, i_c, T2_, Y24_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y24_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_c2, i_c2, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(y,c2,a,c) Y24(c2,w) 
do s_y = 0, nir-1
do s_c2 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c2) == IEOR(s_a,s_c) .and. &
IEOR(s_c2,s_w) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_a, s_c2, s_y)%array(i_a, i_c2, i_y) & 
  * Y24_(s_w, s_c2)%array(i_w, i_c2)
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

end subroutine g_sigma_ccvv_ccvv_no0_x44



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y45(sc1, ic1, V2, Y25, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y45(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y45



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y45(s_c1, i_c1, V2_, Y25_, nir, nsym, psym)

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

integer :: s_c2, i_c2, s_w, i_w
! Y25(c2,w)  <-- 
! (    1.00000000)  V2(c1,c2,c1,w) 
do s_c2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c2,s_w) == 0 .and. & 
IEOR(s_c1,s_c2) == IEOR(s_c1,s_w)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Y25_(s_w, s_c2)%array(i_w, i_c2) =  &
    Y25_(s_w, s_c2)%array(i_w, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_w, s_c1, s_c2)%array(i_w, i_c1, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y45



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x45(sc, ic, T2, Y25, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), Y25(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y25, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x45(sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x45



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x45(s_c, i_c, T2_, Y25_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y25_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_c2, i_c2, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(y,c2,a,c) Y25(c2,w) 
do s_y = 0, nir-1
do s_c2 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c2) == IEOR(s_a,s_c) .and. &
IEOR(s_c2,s_w) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_a, s_c2, s_y)%array(i_a, i_c2, i_y) & 
  * Y25_(s_w, s_c2)%array(i_w, i_c2)
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

end subroutine g_sigma_ccvv_ccvv_no0_x45



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y46(sc1, ic1, V2, Y26, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y46(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y46



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y46(s_c1, i_c1, V2_, Y26_, nir, nsym, psym)

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

integer :: s_c2, i_c2, s_w, i_w
! Y26(c2,w)  <-- 
! (    1.00000000)  V2(c1,c1,c2,w) 
do s_c2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c2,s_w) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c2,s_w)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Y26_(s_w, s_c2)%array(i_w, i_c2) =  &
    Y26_(s_w, s_c2)%array(i_w, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_w, s_c2, s_c1)%array(i_w, i_c2, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y46



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x46(sc, ic, T2, Y26, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), Y26(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y26, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x46(sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x46



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x46(s_c, i_c, T2_, Y26_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y26_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_y, i_y, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -8.00000000) T2(c2,y,a,c) Y26(c2,w) 
do s_c2 = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c2,s_y) == IEOR(s_a,s_c) .and. &
IEOR(s_c2,s_w) == 0) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 8.00000000d+00 & 
  * T2_(s_a, s_y, s_c2)%array(i_a, i_y, i_c2) & 
  * Y26_(s_w, s_c2)%array(i_w, i_c2)
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

end subroutine g_sigma_ccvv_ccvv_no0_x46



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y47(sc1, ic1, V2, Y27, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y27(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y27, nir, nsym, psym) ! -> Ycc (allocate) 
call g_sigma_ccvv_ccvv_y47(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y47



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y47(s_c1, i_c1, V2_, Y27_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y27_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_w, i_w
! Y27(c2,w)  <-- 
! (    1.00000000)  V2(c1,c2,c1,w) 
do s_c2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c2,s_w) == 0 .and. & 
IEOR(s_c1,s_c2) == IEOR(s_c1,s_w)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Y27_(s_w, s_c2)%array(i_w, i_c2) =  &
    Y27_(s_w, s_c2)%array(i_w, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_w, s_c1, s_c2)%array(i_w, i_c1, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y47



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x47(sc, ic, T2, Y27, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), Y27(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y27, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x47(sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x47



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x47(s_c, i_c, T2_, Y27_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y27_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_y, i_y, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(c2,y,a,c) Y27(c2,w) 
do s_c2 = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c2,s_y) == IEOR(s_a,s_c) .and. &
IEOR(s_c2,s_w) == 0) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_a, s_y, s_c2)%array(i_a, i_y, i_c2) & 
  * Y27_(s_w, s_c2)%array(i_w, i_c2)
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

end subroutine g_sigma_ccvv_ccvv_no0_x47



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x48(sc, ic, sc1, ic1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x48(sc, ic, sc1, ic1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x48



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x48(s_c, i_c, s_c1, i_c1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_a, i_a, s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(c2,c1,a,c) V2(c1,y,c2,w) 
do s_c2 = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c2,s_c1) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_y) == IEOR(s_c2,s_w)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_a, s_c1, s_c2)%array(i_a, i_c1, i_c2) & 
  * V2_(s_w, s_c2, s_y)%array(i_w, i_c2, i_y)
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

end subroutine g_sigma_ccvv_ccvv_no0_x48



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x49(sc, ic, sc1, ic1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x49(sc, ic, sc1, ic1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x49



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x49(s_c, i_c, s_c1, i_c1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_a, i_a, s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(c2,c1,a,c) V2(c1,w,c2,y) 
do s_c2 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c2,s_c1) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_w) == IEOR(s_c2,s_y)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_a, s_c1, s_c2)%array(i_a, i_c1, i_c2) & 
  * V2_(s_y, s_c2, s_w)%array(i_y, i_c2, i_w)
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

end subroutine g_sigma_ccvv_ccvv_no0_x49



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y50(sc1, ic1, V2, Y28, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y50(sc1, ic1, h2_i, Y28, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_y50



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y50(s_c1, i_c1, V2_, Y28_, nir, nsym, psym)

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
! (    1.00000000)  V2(c1,c1,c2,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == IEOR(s_c2,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y28_ = Y28_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c2, s_c1)%array(i_c2, i_c2, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y50



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x50(sa, ia, sc, ic, Y28, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: Y28
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x50(sa, ia, sc, ic, Y28, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x50



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x50(s_a, i_a, s_c, i_c, Y28, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y28
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    8.00000000) Y28 T2(y,w,c,a) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_c,s_a)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * Y28 & 
  * T2_(s_c, s_w, s_y)%array(i_c, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x50



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y51(sc1, ic1, V2, Y29, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y29
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_sigma_ccvv_ccvv_y51(sc1, ic1, h2_i, Y29, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_y51



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y51(s_c1, i_c1, V2_, Y29_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y29_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y29()  <-- 
! (    1.00000000)  V2(c1,c2,c1,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c2) == IEOR(s_c1,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y29_ = Y29_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c1, s_c2)%array(i_c2, i_c1, i_c2)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y51



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x51(sa, ia, sc, ic, Y29, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: Y29
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x51(sa, ia, sc, ic, Y29, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x51



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x51(s_a, i_a, s_c, i_c, Y29, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y29
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -4.00000000) Y29 T2(y,w,c,a) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_c,s_a)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * Y29 & 
  * T2_(s_c, s_w, s_y)%array(i_c, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x51



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y52(sc1, ic1, V2, Y30, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y52(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y52



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y52(s_c1, i_c1, V2_, Y30_, nir, nsym, psym)

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

integer :: s_c2, i_c2, s_w, i_w
! Y30(c2,w)  <-- 
! (    1.00000000)  V2(c1,c1,c2,w) 
do s_c2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c2,s_w) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c2,s_w)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Y30_(s_w, s_c2)%array(i_w, i_c2) =  &
    Y30_(s_w, s_c2)%array(i_w, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_w, s_c2, s_c1)%array(i_w, i_c2, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y52



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x52(sa, ia, sc, ic, T2, Y30, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), Y30(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y30, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x52(sa, ia, sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x52



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x52(s_a, i_a, s_c, i_c, T2_, Y30_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: Y30_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_c2, i_c2, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -8.00000000) T2(y,c2,c,a) Y30(c2,w) 
do s_y = 0, nir-1
do s_c2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c2) == IEOR(s_c,s_a) .and. &
IEOR(s_c2,s_w) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 8.00000000d+00 & 
  * T2_(s_c, s_c2, s_y)%array(i_c, i_c2, i_y) & 
  * Y30_(s_w, s_c2)%array(i_w, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x52



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y53(sc1, ic1, V2, Y31, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y53(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y53



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y53(s_c1, i_c1, V2_, Y31_, nir, nsym, psym)

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

integer :: s_c2, i_c2, s_w, i_w
! Y31(c2,w)  <-- 
! (    1.00000000)  V2(c1,c2,c1,w) 
do s_c2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c2,s_w) == 0 .and. & 
IEOR(s_c1,s_c2) == IEOR(s_c1,s_w)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Y31_(s_w, s_c2)%array(i_w, i_c2) =  &
    Y31_(s_w, s_c2)%array(i_w, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_w, s_c1, s_c2)%array(i_w, i_c1, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y53



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x53(sa, ia, sc, ic, T2, Y31, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), Y31(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y31, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x53(sa, ia, sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x53



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x53(s_a, i_a, s_c, i_c, T2_, Y31_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: Y31_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_c2, i_c2, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(y,c2,c,a) Y31(c2,w) 
do s_y = 0, nir-1
do s_c2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c2) == IEOR(s_c,s_a) .and. &
IEOR(s_c2,s_w) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_c, s_c2, s_y)%array(i_c, i_c2, i_y) & 
  * Y31_(s_w, s_c2)%array(i_w, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x53



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y54(sc1, ic1, V2, Y32, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y32
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_sigma_ccvv_ccvv_y54(sc1, ic1, h2_i, Y32, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_y54



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y54(s_c1, i_c1, V2_, Y32_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y32_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y32()  <-- 
! (    1.00000000)  V2(c1,c1,c2,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c1) == IEOR(s_c2,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y32_ = Y32_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c2, s_c1)%array(i_c2, i_c2, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y54



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x54(sa, ia, sc, ic, Y32, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: Y32
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x54(sa, ia, sc, ic, Y32, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x54



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x54(s_a, i_a, s_c, i_c, Y32, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y32
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -4.00000000) Y32 T2(w,y,c,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_c,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * Y32 & 
  * T2_(s_c, s_y, s_w)%array(i_c, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x54



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y55(sc1, ic1, V2, Y33, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y33
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_sigma_ccvv_ccvv_y55(sc1, ic1, h2_i, Y33, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_y55



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y55(s_c1, i_c1, V2_, Y33_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: Y33_
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2
! Y33()  <-- 
! (    1.00000000)  V2(c1,c2,c1,c2) 
do s_c2 = 0, nir-1
if( &
IEOR(s_c1,s_c2) == IEOR(s_c1,s_c2)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
Y33_ = Y33_ &
  + 1.00000000d+00 & 
  * V2_(s_c2, s_c1, s_c2)%array(i_c2, i_c1, i_c2)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y55



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x55(sa, ia, sc, ic, Y33, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: Y33
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x55(sa, ia, sc, ic, Y33, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x55



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x55(s_a, i_a, s_c, i_c, Y33, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y33
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    2.00000000) Y33 T2(w,y,c,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_c,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * Y33 & 
  * T2_(s_c, s_y, s_w)%array(i_c, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x55



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y56(sc1, ic1, V2, Y34, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y34, nir, nsym, psym) ! -> Ycc (allocate) 
call g_sigma_ccvv_ccvv_y56(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y56



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y56(s_c1, i_c1, V2_, Y34_, nir, nsym, psym)

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

integer :: s_c2, i_c2, s_w, i_w
! Y34(c2,w)  <-- 
! (    1.00000000)  V2(c1,c1,c2,w) 
do s_c2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c2,s_w) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c2,s_w)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Y34_(s_w, s_c2)%array(i_w, i_c2) =  &
    Y34_(s_w, s_c2)%array(i_w, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_w, s_c2, s_c1)%array(i_w, i_c2, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y56



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x56(sa, ia, sc, ic, T2, Y34, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), Y34(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y34, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x56(sa, ia, sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x56



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x56(s_a, i_a, s_c, i_c, T2_, Y34_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: Y34_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(c2,y,c,a) Y34(c2,w) 
do s_c2 = 0, nir-1
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c2,s_y) == IEOR(s_c,s_a) .and. &
IEOR(s_c2,s_w) == 0) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_c, s_y, s_c2)%array(i_c, i_y, i_c2) & 
  * Y34_(s_w, s_c2)%array(i_w, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x56



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y57(sc1, ic1, V2, Y35, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y35, nir, nsym, psym) ! -> Ycc (allocate) 
call g_sigma_ccvv_ccvv_y57(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y57



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y57(s_c1, i_c1, V2_, Y35_, nir, nsym, psym)

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

integer :: s_c2, i_c2, s_w, i_w
! Y35(c2,w)  <-- 
! (    1.00000000)  V2(c1,c2,c1,w) 
do s_c2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c2,s_w) == 0 .and. & 
IEOR(s_c1,s_c2) == IEOR(s_c1,s_w)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Y35_(s_w, s_c2)%array(i_w, i_c2) =  &
    Y35_(s_w, s_c2)%array(i_w, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_w, s_c1, s_c2)%array(i_w, i_c1, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y57



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x57(sa, ia, sc, ic, T2, Y35, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), Y35(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y35, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x57(sa, ia, sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x57



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x57(s_a, i_a, s_c, i_c, T2_, Y35_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: Y35_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(c2,y,c,a) Y35(c2,w) 
do s_c2 = 0, nir-1
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c2,s_y) == IEOR(s_c,s_a) .and. &
IEOR(s_c2,s_w) == 0) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_c, s_y, s_c2)%array(i_c, i_y, i_c2) & 
  * Y35_(s_w, s_c2)%array(i_w, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x57



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y58(sc1, ic1, V2, Y36, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y36, nir, nsym, psym) ! -> Ycc (allocate) 
call g_sigma_ccvv_ccvv_y58(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y58



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y58(s_c1, i_c1, V2_, Y36_, nir, nsym, psym)

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

integer :: s_c2, i_c2, s_y, i_y
! Y36(c2,y)  <-- 
! (    1.00000000)  V2(c1,c1,c2,y) 
do s_c2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c2,s_y) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c2,s_y)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y36_(s_y, s_c2)%array(i_y, i_c2) =  &
    Y36_(s_y, s_c2)%array(i_y, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_c2, s_c1)%array(i_y, i_c2, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y58



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x58(sa, ia, sc, ic, T2, Y36, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), Y36(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y36, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x58(sa, ia, sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x58



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x58(s_a, i_a, s_c, i_c, T2_, Y36_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: Y36_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c2, i_c2, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(w,c2,c,a) Y36(c2,y) 
do s_w = 0, nir-1
do s_c2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c2) == IEOR(s_c,s_a) .and. &
IEOR(s_c2,s_y) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_c, s_c2, s_w)%array(i_c, i_c2, i_w) & 
  * Y36_(s_y, s_c2)%array(i_y, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x58



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y59(sc1, ic1, V2, Y37, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y37, nir, nsym, psym) ! -> Ycc (allocate) 
call g_sigma_ccvv_ccvv_y59(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y59



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y59(s_c1, i_c1, V2_, Y37_, nir, nsym, psym)

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

integer :: s_c2, i_c2, s_y, i_y
! Y37(c2,y)  <-- 
! (    1.00000000)  V2(c1,c2,c1,y) 
do s_c2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c2,s_y) == 0 .and. & 
IEOR(s_c1,s_c2) == IEOR(s_c1,s_y)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y37_(s_y, s_c2)%array(i_y, i_c2) =  &
    Y37_(s_y, s_c2)%array(i_y, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_c1, s_c2)%array(i_y, i_c1, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y59



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x59(sa, ia, sc, ic, T2, Y37, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), Y37(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y37, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x59(sa, ia, sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x59



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x59(s_a, i_a, s_c, i_c, T2_, Y37_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: Y37_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c2, i_c2, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(w,c2,c,a) Y37(c2,y) 
do s_w = 0, nir-1
do s_c2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c2) == IEOR(s_c,s_a) .and. &
IEOR(s_c2,s_y) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_c, s_c2, s_w)%array(i_c, i_c2, i_w) & 
  * Y37_(s_y, s_c2)%array(i_y, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x59



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y60(sc1, ic1, V2, Y38, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y38, nir, nsym, psym) ! -> Ycc (allocate) 
call g_sigma_ccvv_ccvv_y60(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y60



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y60(s_c1, i_c1, V2_, Y38_, nir, nsym, psym)

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

integer :: s_c2, i_c2, s_y, i_y
! Y38(c2,y)  <-- 
! (    1.00000000)  V2(c1,c1,c2,y) 
do s_c2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c2,s_y) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c2,s_y)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y38_(s_y, s_c2)%array(i_y, i_c2) =  &
    Y38_(s_y, s_c2)%array(i_y, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_c2, s_c1)%array(i_y, i_c2, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y60



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x60(sa, ia, sc, ic, T2, Y38, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), Y38(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y38, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x60(sa, ia, sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x60



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x60(s_a, i_a, s_c, i_c, T2_, Y38_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: Y38_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -8.00000000) T2(c2,w,c,a) Y38(c2,y) 
do s_c2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c2,s_w) == IEOR(s_c,s_a) .and. &
IEOR(s_c2,s_y) == 0) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 8.00000000d+00 & 
  * T2_(s_c, s_w, s_c2)%array(i_c, i_w, i_c2) & 
  * Y38_(s_y, s_c2)%array(i_y, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x60



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y61(sc1, ic1, V2, Y39, nir, nsym, psym)

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
call set_symblock_Ycc(sleft2, Y39, nir, nsym, psym) ! -> Ycc (allocate) 
call g_sigma_ccvv_ccvv_y61(sc1, ic1, h2_i, Ycc, nir, nsym, psym)

deallocate(h2_i)
deallocate(Ycc)

end subroutine g_if_sigma_ccvv_ccvv_y61



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y61(s_c1, i_c1, V2_, Y39_, nir, nsym, psym)

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

integer :: s_c2, i_c2, s_y, i_y
! Y39(c2,y)  <-- 
! (    1.00000000)  V2(c1,c2,c1,y) 
do s_c2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c2,s_y) == 0 .and. & 
IEOR(s_c1,s_c2) == IEOR(s_c1,s_y)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
Y39_(s_y, s_c2)%array(i_y, i_c2) =  &
    Y39_(s_y, s_c2)%array(i_y, i_c2) &
  + 1.00000000d+00 & 
  * V2_(s_y, s_c1, s_c2)%array(i_y, i_c1, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y61



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x61(sa, ia, sc, ic, T2, Y39, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), Y39(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Ycc(sleft2, Y39, nir, nsym, psym) ! -> Ycc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x61(sa, ia, sc, ic, av2_i, Ycc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Ycc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x61



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x61(s_a, i_a, s_c, i_c, T2_, Y39_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: Y39_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(c2,w,c,a) Y39(c2,y) 
do s_c2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c2,s_w) == IEOR(s_c,s_a) .and. &
IEOR(s_c2,s_y) == 0) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_c, s_w, s_c2)%array(i_c, i_w, i_c2) & 
  * Y39_(s_y, s_c2)%array(i_y, i_c2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x61



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x62(sa, ia, sc, ic, sc1, ic1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x62(sa, ia, sc, ic, sc1, ic1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x62



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x62(s_a, i_a, s_c, i_c, s_c1, i_c1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(c2,c1,c,a) V2(c1,w,c2,y) 
do s_c2 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c2,s_c1) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_w) == IEOR(s_c2,s_y)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_c, s_c1, s_c2)%array(i_c, i_c1, i_c2) & 
  * V2_(s_y, s_c2, s_w)%array(i_y, i_c2, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x62



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x63(sa, ia, sc, ic, sc1, ic1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x63(sa, ia, sc, ic, sc1, ic1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x63



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x63(s_a, i_a, s_c, i_c, s_c1, i_c1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c2, i_c2, s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(c2,c1,c,a) V2(c1,y,c2,w) 
do s_c2 = 0, nir-1
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c2,s_c1) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_y) == IEOR(s_c2,s_w)) then
do i_c2 = psym(I_BEGIN, I_C, s_c2), psym(I_END, I_C, s_c2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_c, s_c1, s_c2)%array(i_c, i_c1, i_c2) & 
  * V2_(s_w, s_c2, s_y)%array(i_w, i_c2, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x63



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y64(sc1, ic1, V2, Y40, nir, nsym, psym)

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
call set_symblock_Yaa(sleft2, Y40, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ccvv_ccvv_y64(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_y64



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y64(s_c1, i_c1, V2_, Y40_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o2, i_o2
! Y40(o1,o2)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y40_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y40_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_o1, s_c1)%array(i_o2, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y64



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x64(Y40, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: Y40(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y40, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = 0
call g_sigma_ccvv_ccvv_no0_x64(Yaa, x, d1, nir, nsym, psym)

deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_no0_x64



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x64(Y40_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y40_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X()  <-- 
! (    1.00000000)  D1(o1,o2) Y40(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y40_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x64



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x64(sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x64(sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x64



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x64(s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    8.00000000) X T2(w,y,a,c) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_a,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * X & 
  * T2_(s_a, s_y, s_w)%array(i_a, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x64



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x65(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x65(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x65



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x65(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_y, i_y
! X(c1,y)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,y,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c1,s_y) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_y) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
X_(s_y)%array(i_y) =  &
    X_(s_y)%array(i_y) &
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

end subroutine g_sigma_ccvv_ccvv_no0_x65



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x65(sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x65(sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x65



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x65(s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(w,c1,a,c) X(c1,y) 
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c1) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_y) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_a, s_c1, s_w)%array(i_a, i_c1, i_w) & 
  * X_(s_y)%array(i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x65



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x66(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x66(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x66



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x66(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_y, i_y
! X(c1,y)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,y,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c1,s_y) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_y) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
X_(s_y)%array(i_y) =  &
    X_(s_y)%array(i_y) &
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

end subroutine g_sigma_ccvv_ccvv_no0_x66



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x66(sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x66(sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x66



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x66(s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(c1,w,a,c) X(c1,y) 
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_w) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_y) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_a, s_w, s_c1)%array(i_a, i_w, i_c1) & 
  * X_(s_y)%array(i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x66



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x67(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x67(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x67



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x67(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_y, i_y
! X(c1,y)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,o1,y,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c1,s_y) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_o1) == IEOR(s_y,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
X_(s_y)%array(i_y) =  &
    X_(s_y)%array(i_y) &
  + 1.00000000d+00 & 
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

end subroutine g_sigma_ccvv_ccvv_no0_x67



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x67(sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x67(sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x67



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x67(s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -1.00000000) T2(c1,w,a,c) X(c1,y) 
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_w) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_y) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 1.00000000d+00 & 
  * T2_(s_a, s_w, s_c1)%array(i_a, i_w, i_c1) & 
  * X_(s_y)%array(i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x67



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x68(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x68(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x68



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x68(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! X(c1,w)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,w,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c1,s_w) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_w) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
X_(s_w)%array(i_w) =  &
    X_(s_w)%array(i_w) &
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

end subroutine g_sigma_ccvv_ccvv_no0_x68



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x68(sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x68(sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x68



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x68(s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(y,c1,a,c) X(c1,w) 
do s_y = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c1) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_w) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_a, s_c1, s_y)%array(i_a, i_c1, i_y) & 
  * X_(s_w)%array(i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x68



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x69(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x69(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x69



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x69(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! X(c1,w)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,o1,w,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c1,s_w) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_o1) == IEOR(s_w,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
X_(s_w)%array(i_w) =  &
    X_(s_w)%array(i_w) &
  + 1.00000000d+00 & 
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

end subroutine g_sigma_ccvv_ccvv_no0_x69



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x69(sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x69(sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x69



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x69(s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -1.00000000) T2(y,c1,a,c) X(c1,w) 
do s_y = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c1) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_w) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 1.00000000d+00 & 
  * T2_(s_a, s_c1, s_y)%array(i_a, i_c1, i_y) & 
  * X_(s_w)%array(i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x69



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x70(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x70(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x70



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x70(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! X(c1,w)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,w,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c1,s_w) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_w) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
X_(s_w)%array(i_w) =  &
    X_(s_w)%array(i_w) &
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

end subroutine g_sigma_ccvv_ccvv_no0_x70



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x70(sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x70(sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x70



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x70(s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(c1,y,a,c) X(c1,w) 
do s_y = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_y) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_w) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_a, s_y, s_c1)%array(i_a, i_y, i_c1) & 
  * X_(s_w)%array(i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x70



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y71(sc1, ic1, V2, Y41, nir, nsym, psym)

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
call set_symblock_Yaa(sleft2, Y41, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ccvv_ccvv_y71(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_y71



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y71(s_c1, i_c1, V2_, Y41_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o2, i_o2
! Y41(o1,o2)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y41_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y41_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_o1, s_c1)%array(i_o2, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y71



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x71(Y41, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: Y41(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y41, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = 0
call g_sigma_ccvv_ccvv_no0_x71(Yaa, x, d1, nir, nsym, psym)

deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_no0_x71



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x71(Y41_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y41_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X()  <-- 
! (    1.00000000)  D1(o1,o2) Y41(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y41_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x71



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x71(sa, ia, sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x71(sa, ia, sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x71



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x71(s_a, i_a, s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    8.00000000) X T2(y,w,c,a) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_c,s_a)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * X & 
  * T2_(s_c, s_w, s_y)%array(i_c, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x71



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x72(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x72(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x72



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x72(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! X(c1,w)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,w,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c1,s_w) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_w) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
X_(s_w)%array(i_w) =  &
    X_(s_w)%array(i_w) &
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

end subroutine g_sigma_ccvv_ccvv_no0_x72



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x72(sa, ia, sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x72(sa, ia, sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x72



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x72(s_a, i_a, s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(y,c1,c,a) X(c1,w) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c1) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_w) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_c, s_c1, s_y)%array(i_c, i_c1, i_y) & 
  * X_(s_w)%array(i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x72



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x73(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x73(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x73



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x73(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! X(c1,w)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,w,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c1,s_w) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_w) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
X_(s_w)%array(i_w) =  &
    X_(s_w)%array(i_w) &
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

end subroutine g_sigma_ccvv_ccvv_no0_x73



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x73(sa, ia, sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x73(sa, ia, sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x73



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x73(s_a, i_a, s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(c1,y,c,a) X(c1,w) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_y) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_w) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_c, s_y, s_c1)%array(i_c, i_y, i_c1) & 
  * X_(s_w)%array(i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x73



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x74(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x74(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x74



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x74(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! X(c1,w)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,o1,w,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c1,s_w) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_o1) == IEOR(s_w,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
X_(s_w)%array(i_w) =  &
    X_(s_w)%array(i_w) &
  + 1.00000000d+00 & 
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

end subroutine g_sigma_ccvv_ccvv_no0_x74



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x74(sa, ia, sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x74(sa, ia, sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x74



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x74(s_a, i_a, s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -1.00000000) T2(c1,y,c,a) X(c1,w) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_y) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_w) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 1.00000000d+00 & 
  * T2_(s_c, s_y, s_c1)%array(i_c, i_y, i_c1) & 
  * X_(s_w)%array(i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x74



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x75(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x75(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x75



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x75(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_y, i_y
! X(c1,y)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,y,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c1,s_y) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_y) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
X_(s_y)%array(i_y) =  &
    X_(s_y)%array(i_y) &
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

end subroutine g_sigma_ccvv_ccvv_no0_x75



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x75(sa, ia, sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x75(sa, ia, sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x75



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x75(s_a, i_a, s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(w,c1,c,a) X(c1,y) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c1) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_y) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_c, s_c1, s_w)%array(i_c, i_c1, i_w) & 
  * X_(s_y)%array(i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x75



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x76(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x76(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x76



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x76(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_y, i_y
! X(c1,y)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,o1,y,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c1,s_y) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_o1) == IEOR(s_y,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
X_(s_y)%array(i_y) =  &
    X_(s_y)%array(i_y) &
  + 1.00000000d+00 & 
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

end subroutine g_sigma_ccvv_ccvv_no0_x76



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x76(sa, ia, sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x76(sa, ia, sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x76



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x76(s_a, i_a, s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -1.00000000) T2(w,c1,c,a) X(c1,y) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c1) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_y) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 1.00000000d+00 & 
  * T2_(s_c, s_c1, s_w)%array(i_c, i_c1, i_w) & 
  * X_(s_y)%array(i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x76



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x77(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x77(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x77



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x77(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_y, i_y
! X(c1,y)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,y,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c1,s_y) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_y) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
X_(s_y)%array(i_y) =  &
    X_(s_y)%array(i_y) &
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

end subroutine g_sigma_ccvv_ccvv_no0_x77



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x77(sa, ia, sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x77(sa, ia, sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x77



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x77(s_a, i_a, s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(c1,w,c,a) X(c1,y) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_w) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_y) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_c, s_w, s_c1)%array(i_c, i_w, i_c1) & 
  * X_(s_y)%array(i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x77



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y78(sc1, ic1, V2, Y42, nir, nsym, psym)

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
call set_symblock_Yaa(sleft2, Y42, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ccvv_ccvv_y78(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_y78



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y78(s_c1, i_c1, V2_, Y42_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o2, i_o2
! Y42(o1,o2)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y42_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y42_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_c1, s_o1)%array(i_o2, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y78



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x78(Y42, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: Y42(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y42, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = 0
call g_sigma_ccvv_ccvv_no0_x78(Yaa, x, d1, nir, nsym, psym)

deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_no0_x78



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x78(Y42_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y42_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X()  <-- 
! (    1.00000000)  D1(o1,o2) Y42(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y42_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x78



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x78(sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x78(sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x78



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x78(s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) X T2(w,y,a,c) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_a,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * X & 
  * T2_(s_a, s_y, s_w)%array(i_a, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x78



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x79(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x79(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x79



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x79(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_y, i_y
! X(c1,y)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,o1,y,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c1,s_y) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_o1) == IEOR(s_y,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
X_(s_y)%array(i_y) =  &
    X_(s_y)%array(i_y) &
  + 1.00000000d+00 & 
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

end subroutine g_sigma_ccvv_ccvv_no0_x79



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x79(sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x79(sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x79



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x79(s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(w,c1,a,c) X(c1,y) 
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c1) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_y) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_a, s_c1, s_w)%array(i_a, i_c1, i_w) & 
  * X_(s_y)%array(i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x79



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x80(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x80(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x80



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x80(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! X(c1,w)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,o1,w,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c1,s_w) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_o1) == IEOR(s_w,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
X_(s_w)%array(i_w) =  &
    X_(s_w)%array(i_w) &
  + 1.00000000d+00 & 
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

end subroutine g_sigma_ccvv_ccvv_no0_x80



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x80(sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x80(sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x80



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x80(s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(c1,y,a,c) X(c1,w) 
do s_y = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_y) == IEOR(s_a,s_c) .and. &
IEOR(s_c1,s_w) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_a, s_y, s_c1)%array(i_a, i_y, i_c1) & 
  * X_(s_w)%array(i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x80



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y81(sc1, ic1, V2, Y43, nir, nsym, psym)

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
call set_symblock_Yaa(sleft2, Y43, nir, nsym, psym) ! -> Yaa (allocate) 
call g_sigma_ccvv_ccvv_y81(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_y81



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y81(s_c1, i_c1, V2_, Y43_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o2, i_o2
! Y43(o1,o2)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Y43_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    Y43_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o2, s_c1, s_o1)%array(i_o2, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y81



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x81(Y43, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: Y43(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y43, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = 0
call g_sigma_ccvv_ccvv_no0_x81(Yaa, x, d1, nir, nsym, psym)

deallocate(Yaa)

end subroutine g_if_sigma_ccvv_ccvv_no0_x81



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x81(Y43_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y43_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X()  <-- 
! (    1.00000000)  D1(o1,o2) Y43(o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_o1,s_o2) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * Y43_(s_o2, s_o1)%array(i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x81



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x81(sa, ia, sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x81(sa, ia, sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x81



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x81(s_a, i_a, s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -4.00000000) X T2(y,w,c,a) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_c,s_a)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * X & 
  * T2_(s_c, s_w, s_y)%array(i_c, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x81



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x82(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x82(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x82



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x82(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_w, i_w
! X(c1,w)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,o1,w,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_c1,s_w) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_o1) == IEOR(s_w,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
X_(s_w)%array(i_w) =  &
    X_(s_w)%array(i_w) &
  + 1.00000000d+00 & 
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

end subroutine g_sigma_ccvv_ccvv_no0_x82



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x82(sa, ia, sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x82(sa, ia, sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x82



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x82(s_a, i_a, s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(y,c1,c,a) X(c1,w) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c1) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_w) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_c, s_c1, s_y)%array(i_c, i_c1, i_y) & 
  * X_(s_w)%array(i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x82



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x83(sc1, ic1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call g_sigma_ccvv_ccvv_no0_x83(sc1, ic1, h2_i, xc, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xc)

end subroutine g_if_sigma_ccvv_ccvv_no0_x83



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x83(s_c1, i_c1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_y, i_y
! X(c1,y)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c1,o1,y,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_c1,s_y) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c1,s_o1) == IEOR(s_y,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
X_(s_y)%array(i_y) =  &
    X_(s_y)%array(i_y) &
  + 1.00000000d+00 & 
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

end subroutine g_sigma_ccvv_ccvv_no0_x83



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x83(sa, ia, sc, ic, sc1, ic1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic, sc1, ic1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc1

call set_symblock_xc(sleft, x, nir, nsym, psym) ! -> xc (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x83(sa, ia, sc, ic, sc1, ic1, av2_i, xc, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xc)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x83



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x83(s_a, i_a, s_c, i_c, s_c1, i_c1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(c1,w,c,a) X(c1,y) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_w) == IEOR(s_c,s_a) .and. &
IEOR(s_c1,s_y) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_c, s_w, s_c1)%array(i_c, i_w, i_c1) & 
  * X_(s_y)%array(i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x83



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y84(sc1, ic1, V2, Y44, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y44, nir, nsym, psym) ! -> Yvv (allocate) 
call g_sigma_ccvv_ccvv_y84(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y84



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y84(s_c1, i_c1, V2_, Y44_, nir, nsym, psym)

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

integer :: s_c, i_c, s_v1, i_v1
! Y44(c,v1)  <-- 
! (    1.00000000)  V2(c1,c1,c,v1) 
do s_c = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c,s_v1)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y44_(s_v1, s_c)%array(i_v1, i_c) =  &
    Y44_(s_v1, s_c)%array(i_v1, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c, s_c1)%array(i_v1, i_c, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y84



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x84(sa, ia, sc, ic, T2, Y44, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), Y44(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y44, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x84(sa, ia, sc, ic, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x84



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x84(s_a, i_a, s_c, i_c, T2_, Y44_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: Y44_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_v1, i_v1
! S2(w,y,a,c)  <-- 
! (    8.00000000) T2(y,w,v1,a) Y44(c,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * T2_(s_v1, s_w, s_y)%array(i_v1, i_w, i_y) & 
  * Y44_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x84



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x85(sa, ia, sc, ic, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x85(sa, ia, sc, ic, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x85



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x85(s_a, i_a, s_c, i_c, T2_, V2_, S2_, nir, nsym, psym)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_c1, i_c1, s_v1, i_v1, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(y,c1,v1,a) V2(c,v1,c1,w) 
do s_y = 0, nir-1
do s_c1 = 0, nir-1
do s_v1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c1) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == IEOR(s_c1,s_w)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_v1, s_c1, s_y)%array(i_v1, i_c1, i_y) & 
  * V2_(s_w, s_c1, s_v1)%array(i_w, i_c1, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x85



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x86(sa, ia, sc, ic, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x86(sa, ia, sc, ic, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x86



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x86(s_a, i_a, s_c, i_c, T2_, V2_, S2_, nir, nsym, psym)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_y, i_y, s_v1, i_v1, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(c1,y,v1,a) V2(c,v1,c1,w) 
do s_c1 = 0, nir-1
do s_y = 0, nir-1
do s_v1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_y) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == IEOR(s_c1,s_w)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_v1, s_y, s_c1)%array(i_v1, i_y, i_c1) & 
  * V2_(s_w, s_c1, s_v1)%array(i_w, i_c1, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x86



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x87(sa, ia, sc, ic, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x87(sa, ia, sc, ic, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x87



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x87(s_a, i_a, s_c, i_c, T2_, V2_, S2_, nir, nsym, psym)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_y, i_y, s_v1, i_v1, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(c1,y,v1,a) V2(c,w,c1,v1) 
do s_c1 = 0, nir-1
do s_y = 0, nir-1
do s_v1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_y) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_w) == IEOR(s_c1,s_v1)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_v1, s_y, s_c1)%array(i_v1, i_y, i_c1) & 
  * V2_(s_v1, s_c1, s_w)%array(i_v1, i_c1, i_w)
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

end subroutine g_sigma_ccvv_ccvv_no0_x87



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x88(sa, ia, sc, ic, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x88(sa, ia, sc, ic, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x88



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x88(s_a, i_a, s_c, i_c, T2_, V2_, S2_, nir, nsym, psym)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c1, i_c1, s_v1, i_v1, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(w,c1,v1,a) V2(c,v1,c1,y) 
do s_w = 0, nir-1
do s_c1 = 0, nir-1
do s_v1 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c1) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == IEOR(s_c1,s_y)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_v1, s_c1, s_w)%array(i_v1, i_c1, i_w) & 
  * V2_(s_y, s_c1, s_v1)%array(i_y, i_c1, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x88



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x89(sa, ia, sc, ic, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x89(sa, ia, sc, ic, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x89



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x89(s_a, i_a, s_c, i_c, T2_, V2_, S2_, nir, nsym, psym)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c1, i_c1, s_v1, i_v1, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(w,c1,v1,a) V2(c,y,c1,v1) 
do s_w = 0, nir-1
do s_c1 = 0, nir-1
do s_v1 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c1) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_y) == IEOR(s_c1,s_v1)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_v1, s_c1, s_w)%array(i_v1, i_c1, i_w) & 
  * V2_(s_v1, s_c1, s_y)%array(i_v1, i_c1, i_y)
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

end subroutine g_sigma_ccvv_ccvv_no0_x89



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x90(sa, ia, sc, ic, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x90(sa, ia, sc, ic, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x90



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x90(s_a, i_a, s_c, i_c, T2_, V2_, S2_, nir, nsym, psym)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_w, i_w, s_v1, i_v1, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(c1,w,v1,a) V2(c,v1,c1,y) 
do s_c1 = 0, nir-1
do s_w = 0, nir-1
do s_v1 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_w) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == IEOR(s_c1,s_y)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_v1, s_w, s_c1)%array(i_v1, i_w, i_c1) & 
  * V2_(s_y, s_c1, s_v1)%array(i_y, i_c1, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x90



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y91(sc1, ic1, V2, Y45, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y45, nir, nsym, psym) ! -> Yvv (allocate) 
call g_sigma_ccvv_ccvv_y91(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y91



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y91(s_c1, i_c1, V2_, Y45_, nir, nsym, psym)

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

integer :: s_a, i_a, s_v1, i_v1
! Y45(a,v1)  <-- 
! (    1.00000000)  V2(c1,c1,a,v1) 
do s_a = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_v1)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y45_(s_v1, s_a)%array(i_v1, i_a) =  &
    Y45_(s_v1, s_a)%array(i_v1, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_a, s_c1)%array(i_v1, i_a, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y91



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x91(sc, ic, T2, Y45, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), Y45(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y45, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x91(sc, ic, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x91



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x91(s_c, i_c, T2_, Y45_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y45_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_v1, i_v1, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    8.00000000) T2(w,y,v1,c) Y45(a,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_v1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * T2_(s_v1, s_y, s_w)%array(i_v1, i_y, i_w) & 
  * Y45_(s_v1, s_a)%array(i_v1, i_a)
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

end subroutine g_sigma_ccvv_ccvv_no0_x91



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x92(sa, ia, sc, ic, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x92(sa, ia, sc, ic, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x92



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x92(s_a, i_a, s_c, i_c, T2_, V2_, S2_, nir, nsym, psym)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c1, i_c1, s_v1, i_v1, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(w,c1,v1,c) V2(a,v1,c1,y) 
do s_w = 0, nir-1
do s_c1 = 0, nir-1
do s_v1 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c1) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == IEOR(s_c1,s_y)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_v1, s_c1, s_w)%array(i_v1, i_c1, i_w) & 
  * V2_(s_y, s_c1, s_v1)%array(i_y, i_c1, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x92



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x93(sa, ia, sc, ic, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x93(sa, ia, sc, ic, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x93



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x93(s_a, i_a, s_c, i_c, T2_, V2_, S2_, nir, nsym, psym)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_w, i_w, s_v1, i_v1, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(c1,w,v1,c) V2(a,v1,c1,y) 
do s_c1 = 0, nir-1
do s_w = 0, nir-1
do s_v1 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_w) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == IEOR(s_c1,s_y)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_v1, s_w, s_c1)%array(i_v1, i_w, i_c1) & 
  * V2_(s_y, s_c1, s_v1)%array(i_y, i_c1, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x93



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x94(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x94(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x94



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x94(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(c1,w,v1,c) V2(v1,c1,y,a) 
do s_c1 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_w) == IEOR(s_v1,s_c) .and. &
IEOR(s_v1,s_c1) == IEOR(s_y,s_a)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_v1, s_w, s_c1)%array(i_v1, i_w, i_c1) & 
  * V2_(s_a, s_y, s_c1)%array(i_a, i_y, i_c1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x94



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x95(sa, ia, sc, ic, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x95(sa, ia, sc, ic, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x95



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x95(s_a, i_a, s_c, i_c, T2_, V2_, S2_, nir, nsym, psym)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_c1, i_c1, s_v1, i_v1, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(y,c1,v1,c) V2(a,v1,c1,w) 
do s_y = 0, nir-1
do s_c1 = 0, nir-1
do s_v1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c1) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == IEOR(s_c1,s_w)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_v1, s_c1, s_y)%array(i_v1, i_c1, i_y) & 
  * V2_(s_w, s_c1, s_v1)%array(i_w, i_c1, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x95



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x96(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x96(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x96



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x96(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_c1, i_c1, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(y,c1,v1,c) V2(v1,c1,w,a) 
do s_y = 0, nir-1
do s_c1 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c1) == IEOR(s_v1,s_c) .and. &
IEOR(s_v1,s_c1) == IEOR(s_w,s_a)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_v1, s_c1, s_y)%array(i_v1, i_c1, i_y) & 
  * V2_(s_a, s_w, s_c1)%array(i_a, i_w, i_c1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x96



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x97(sa, ia, sc, ic, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x97(sa, ia, sc, ic, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x97



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x97(s_a, i_a, s_c, i_c, T2_, V2_, S2_, nir, nsym, psym)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_y, i_y, s_v1, i_v1, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(c1,y,v1,c) V2(a,v1,c1,w) 
do s_c1 = 0, nir-1
do s_y = 0, nir-1
do s_v1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_y) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == IEOR(s_c1,s_w)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_v1, s_y, s_c1)%array(i_v1, i_y, i_c1) & 
  * V2_(s_w, s_c1, s_v1)%array(i_w, i_c1, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x97



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y98(sc1, ic1, V2, Y46, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y46, nir, nsym, psym) ! -> Yvv (allocate) 
call g_sigma_ccvv_ccvv_y98(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y98



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y98(s_c1, i_c1, V2_, Y46_, nir, nsym, psym)

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

integer :: s_c, i_c, s_v1, i_v1
! Y46(c,v1)  <-- 
! (    1.00000000)  V2(c1,c1,c,v1) 
do s_c = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c,s_v1)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y46_(s_v1, s_c)%array(i_v1, i_c) =  &
    Y46_(s_v1, s_c)%array(i_v1, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c, s_c1)%array(i_v1, i_c, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y98



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x98(sc, ic, sv1, iv1, T2, Y46, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y46(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y46, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x98(sc, ic, sv1, iv1, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x98



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x98(s_c, i_c, s_v1, i_v1, T2_, Y46_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y46_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    8.00000000) T2(w,y,a,v1) Y46(c,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * T2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) & 
  * Y46_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x98



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x99(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x99(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x99



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x99(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c1, i_c1, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(w,c1,a,v1) V2(c,v1,c1,y) 
do s_w = 0, nir-1
do s_c1 = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c1) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == IEOR(s_c1,s_y)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_a, s_c1, s_w)%array(i_a, i_c1, i_w) & 
  * V2_(s_y, s_c1, s_v1)%array(i_y, i_c1, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x99



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x100(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x100(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x100



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x100(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_w, i_w, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(c1,w,a,v1) V2(c,v1,c1,y) 
do s_c1 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_w) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == IEOR(s_c1,s_y)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_a, s_w, s_c1)%array(i_a, i_w, i_c1) & 
  * V2_(s_y, s_c1, s_v1)%array(i_y, i_c1, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x100



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x101(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x101(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x101



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x101(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_w, i_w, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(c1,w,a,v1) V2(c,y,c1,v1) 
do s_c1 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_w) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_y) == IEOR(s_c1,s_v1)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_a, s_w, s_c1)%array(i_a, i_w, i_c1) & 
  * V2_(s_v1, s_c1, s_y)%array(i_v1, i_c1, i_y)
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

end subroutine g_sigma_ccvv_ccvv_no0_x101



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x102(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x102(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x102



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x102(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_c1, i_c1, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(y,c1,a,v1) V2(c,v1,c1,w) 
do s_y = 0, nir-1
do s_c1 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c1) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == IEOR(s_c1,s_w)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_a, s_c1, s_y)%array(i_a, i_c1, i_y) & 
  * V2_(s_w, s_c1, s_v1)%array(i_w, i_c1, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x102



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x103(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x103(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x103



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x103(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_c1, i_c1, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(y,c1,a,v1) V2(c,w,c1,v1) 
do s_y = 0, nir-1
do s_c1 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c1) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_w) == IEOR(s_c1,s_v1)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_a, s_c1, s_y)%array(i_a, i_c1, i_y) & 
  * V2_(s_v1, s_c1, s_w)%array(i_v1, i_c1, i_w)
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

end subroutine g_sigma_ccvv_ccvv_no0_x103



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x104(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x104(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x104



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x104(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_y, i_y, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(c1,y,a,v1) V2(c,v1,c1,w) 
do s_c1 = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_y) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == IEOR(s_c1,s_w)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_a, s_y, s_c1)%array(i_a, i_y, i_c1) & 
  * V2_(s_w, s_c1, s_v1)%array(i_w, i_c1, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x104



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y105(sc1, ic1, V2, Y47, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y47, nir, nsym, psym) ! -> Yvv (allocate) 
call g_sigma_ccvv_ccvv_y105(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y105



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y105(s_c1, i_c1, V2_, Y47_, nir, nsym, psym)

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

integer :: s_a, i_a, s_v1, i_v1
! Y47(a,v1)  <-- 
! (    1.00000000)  V2(c1,c1,a,v1) 
do s_a = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_v1)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y47_(s_v1, s_a)%array(i_v1, i_a) =  &
    Y47_(s_v1, s_a)%array(i_v1, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_a, s_c1)%array(i_v1, i_a, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y105



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x105(sc, ic, sv1, iv1, T2, Y47, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y47(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y47, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x105(sc, ic, sv1, iv1, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x105



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x105(s_c, i_c, s_v1, i_v1, T2_, Y47_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y47_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    8.00000000) T2(y,w,c,v1) Y47(a,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_c,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * T2_(s_c, s_w, s_y)%array(i_c, i_w, i_y) & 
  * Y47_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x105



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x106(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x106(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x106



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x106(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_c1, i_c1, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(y,c1,c,v1) V2(v1,a,c1,w) 
do s_y = 0, nir-1
do s_c1 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c1) == IEOR(s_c,s_v1) .and. &
IEOR(s_v1,s_a) == IEOR(s_c1,s_w)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_c, s_c1, s_y)%array(i_c, i_c1, i_y) & 
  * V2_(s_w, s_c1, s_a)%array(i_w, i_c1, i_a)
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

end subroutine g_sigma_ccvv_ccvv_no0_x106



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x107(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x107(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x107



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x107(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_y, i_y, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(c1,y,c,v1) V2(v1,a,c1,w) 
do s_c1 = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_y) == IEOR(s_c,s_v1) .and. &
IEOR(s_v1,s_a) == IEOR(s_c1,s_w)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_c, s_y, s_c1)%array(i_c, i_y, i_c1) & 
  * V2_(s_w, s_c1, s_a)%array(i_w, i_c1, i_a)
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

end subroutine g_sigma_ccvv_ccvv_no0_x107



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x108(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x108(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x108



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x108(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(c1,y,c,v1) V2(v1,c1,w,a) 
do s_c1 = 0, nir-1
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_y) == IEOR(s_c,s_v1) .and. &
IEOR(s_v1,s_c1) == IEOR(s_w,s_a)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_c, s_y, s_c1)%array(i_c, i_y, i_c1) & 
  * V2_(s_a, s_w, s_c1)%array(i_a, i_w, i_c1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x108



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x109(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x109(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x109



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x109(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c1, i_c1, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(w,c1,c,v1) V2(v1,a,c1,y) 
do s_w = 0, nir-1
do s_c1 = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c1) == IEOR(s_c,s_v1) .and. &
IEOR(s_v1,s_a) == IEOR(s_c1,s_y)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_c, s_c1, s_w)%array(i_c, i_c1, i_w) & 
  * V2_(s_y, s_c1, s_a)%array(i_y, i_c1, i_a)
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

end subroutine g_sigma_ccvv_ccvv_no0_x109



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x110(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x110(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x110



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x110(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c1, i_c1, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(w,c1,c,v1) V2(v1,c1,y,a) 
do s_w = 0, nir-1
do s_c1 = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c1) == IEOR(s_c,s_v1) .and. &
IEOR(s_v1,s_c1) == IEOR(s_y,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_c, s_c1, s_w)%array(i_c, i_c1, i_w) & 
  * V2_(s_a, s_y, s_c1)%array(i_a, i_y, i_c1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x110



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x111(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x111(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x111



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x111(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_w, i_w, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(c1,w,c,v1) V2(v1,a,c1,y) 
do s_c1 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_w) == IEOR(s_c,s_v1) .and. &
IEOR(s_v1,s_a) == IEOR(s_c1,s_y)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_c, s_w, s_c1)%array(i_c, i_w, i_c1) & 
  * V2_(s_y, s_c1, s_a)%array(i_y, i_c1, i_a)
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

end subroutine g_sigma_ccvv_ccvv_no0_x111



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x112(sa, ia, sc, ic, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x112(sa, ia, sc, ic, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x112



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x112(s_a, i_a, s_c, i_c, T2_, V2_, S2_, nir, nsym, psym)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_w, i_w, s_v1, i_v1, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    8.00000000) T2(c1,w,v1,a) V2(c,y,c1,v1) 
do s_c1 = 0, nir-1
do s_w = 0, nir-1
do s_v1 = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_w) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_y) == IEOR(s_c1,s_v1)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * T2_(s_v1, s_w, s_c1)%array(i_v1, i_w, i_c1) & 
  * V2_(s_v1, s_c1, s_y)%array(i_v1, i_c1, i_y)
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

end subroutine g_sigma_ccvv_ccvv_no0_x112



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x113(sa, ia, sc, ic, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x113(sa, ia, sc, ic, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x113



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x113(s_a, i_a, s_c, i_c, T2_, V2_, S2_, nir, nsym, psym)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_c1, i_c1, s_v1, i_v1, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(y,c1,v1,a) V2(c,w,c1,v1) 
do s_y = 0, nir-1
do s_c1 = 0, nir-1
do s_v1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c1) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_w) == IEOR(s_c1,s_v1)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_v1, s_c1, s_y)%array(i_v1, i_c1, i_y) & 
  * V2_(s_v1, s_c1, s_w)%array(i_v1, i_c1, i_w)
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

end subroutine g_sigma_ccvv_ccvv_no0_x113



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y114(sc1, ic1, V2, Y48, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y114(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y114



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y114(s_c1, i_c1, V2_, Y48_, nir, nsym, psym)

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

integer :: s_c, i_c, s_v1, i_v1
! Y48(c,v1)  <-- 
! (    1.00000000)  V2(c1,c,c1,v1) 
do s_c = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_c1,s_c) == IEOR(s_c1,s_v1)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y48_(s_v1, s_c)%array(i_v1, i_c) =  &
    Y48_(s_v1, s_c)%array(i_v1, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_c)%array(i_v1, i_c1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y114



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x114(sa, ia, sc, ic, T2, Y48, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), Y48(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y48, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x114(sa, ia, sc, ic, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x114



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x114(s_a, i_a, s_c, i_c, T2_, Y48_, S2_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: Y48_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_v1, i_v1
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(y,w,v1,a) Y48(c,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_v1, s_w, s_y)%array(i_v1, i_w, i_y) & 
  * Y48_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x114



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x115(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x115(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x115



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x115(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    8.00000000) T2(c1,y,v1,c) V2(v1,c1,w,a) 
do s_c1 = 0, nir-1
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_y) == IEOR(s_v1,s_c) .and. &
IEOR(s_v1,s_c1) == IEOR(s_w,s_a)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * T2_(s_v1, s_y, s_c1)%array(i_v1, i_y, i_c1) & 
  * V2_(s_a, s_w, s_c1)%array(i_a, i_w, i_c1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x115



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x116(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x116(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x116



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x116(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c1, i_c1, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(w,c1,v1,c) V2(v1,c1,y,a) 
do s_w = 0, nir-1
do s_c1 = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c1) == IEOR(s_v1,s_c) .and. &
IEOR(s_v1,s_c1) == IEOR(s_y,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_v1, s_c1, s_w)%array(i_v1, i_c1, i_w) & 
  * V2_(s_a, s_y, s_c1)%array(i_a, i_y, i_c1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x116



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y117(sc1, ic1, V2, Y49, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y117(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y117



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y117(s_c1, i_c1, V2_, Y49_, nir, nsym, psym)

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

integer :: s_a, i_a, s_v1, i_v1
! Y49(a,v1)  <-- 
! (    1.00000000)  V2(c1,a,c1,v1) 
do s_a = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_v1)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y49_(s_v1, s_a)%array(i_v1, i_a) =  &
    Y49_(s_v1, s_a)%array(i_v1, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_a)%array(i_v1, i_c1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y117



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x117(sc, ic, T2, Y49, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), Y49(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y49, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x117(sc, ic, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x117



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x117(s_c, i_c, T2_, Y49_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y49_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_v1, i_v1, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(w,y,v1,c) Y49(a,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_v1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_v1, s_y, s_w)%array(i_v1, i_y, i_w) & 
  * Y49_(s_v1, s_a)%array(i_v1, i_a)
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

end subroutine g_sigma_ccvv_ccvv_no0_x117



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x118(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x118(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x118



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x118(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c1, i_c1, s_a, i_a, s_y, i_y
! S2(w,y,a,c)  <-- 
! (    8.00000000) T2(w,c1,a,v1) V2(c,y,c1,v1) 
do s_w = 0, nir-1
do s_c1 = 0, nir-1
do s_a = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_c1) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_y) == IEOR(s_c1,s_v1)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * T2_(s_a, s_c1, s_w)%array(i_a, i_c1, i_w) & 
  * V2_(s_v1, s_c1, s_y)%array(i_v1, i_c1, i_y)
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

end subroutine g_sigma_ccvv_ccvv_no0_x118



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y119(sc1, ic1, V2, Y50, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y119(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y119



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y119(s_c1, i_c1, V2_, Y50_, nir, nsym, psym)

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

integer :: s_c, i_c, s_v1, i_v1
! Y50(c,v1)  <-- 
! (    1.00000000)  V2(c1,c,c1,v1) 
do s_c = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_c1,s_c) == IEOR(s_c1,s_v1)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y50_(s_v1, s_c)%array(i_v1, i_c) =  &
    Y50_(s_v1, s_c)%array(i_v1, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_c)%array(i_v1, i_c1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y119



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x119(sc, ic, sv1, iv1, T2, Y50, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y50(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y50, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x119(sc, ic, sv1, iv1, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x119



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x119(s_c, i_c, s_v1, i_v1, T2_, Y50_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y50_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(w,y,a,v1) Y50(c,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) & 
  * Y50_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x119



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x120(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x120(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x120



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x120(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_y, i_y, s_a, i_a, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(c1,y,a,v1) V2(c,w,c1,v1) 
do s_c1 = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_y) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_w) == IEOR(s_c1,s_v1)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_a, s_y, s_c1)%array(i_a, i_y, i_c1) & 
  * V2_(s_v1, s_c1, s_w)%array(i_v1, i_c1, i_w)
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

end subroutine g_sigma_ccvv_ccvv_no0_x120



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x121(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x121(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x121



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x121(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_c1, i_c1, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    8.00000000) T2(y,c1,c,v1) V2(v1,c1,w,a) 
do s_y = 0, nir-1
do s_c1 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_c1) == IEOR(s_c,s_v1) .and. &
IEOR(s_v1,s_c1) == IEOR(s_w,s_a)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 8.00000000d+00 & 
  * T2_(s_c, s_c1, s_y)%array(i_c, i_c1, i_y) & 
  * V2_(s_a, s_w, s_c1)%array(i_a, i_w, i_c1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x121



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_y122(sc1, ic1, V2, Y51, nir, nsym, psym)

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
call g_sigma_ccvv_ccvv_y122(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_ccvv_ccvv_y122



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_y122(s_c1, i_c1, V2_, Y51_, nir, nsym, psym)

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

integer :: s_a, i_a, s_v1, i_v1
! Y51(a,v1)  <-- 
! (    1.00000000)  V2(c1,a,c1,v1) 
do s_a = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_v1)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y51_(s_v1, s_a)%array(i_v1, i_a) =  &
    Y51_(s_v1, s_a)%array(i_v1, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_a)%array(i_v1, i_c1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_y122



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x122(sc, ic, sv1, iv1, T2, Y51, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y51(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y51, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x122(sc, ic, sv1, iv1, av2_i, Yvv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x122



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x122(s_c, i_c, s_v1, i_v1, T2_, Y51_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y51_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -4.00000000) T2(y,w,c,v1) Y51(a,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_c,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 4.00000000d+00 & 
  * T2_(s_c, s_w, s_y)%array(i_c, i_w, i_y) & 
  * Y51_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x122



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x123(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x123(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x123



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x123(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    2.00000000) T2(c1,w,c,v1) V2(v1,c1,y,a) 
do s_c1 = 0, nir-1
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c1,s_w) == IEOR(s_c,s_v1) .and. &
IEOR(s_v1,s_c1) == IEOR(s_y,s_a)) then
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * T2_(s_c, s_w, s_c1)%array(i_c, i_w, i_c1) & 
  * V2_(s_a, s_y, s_c1)%array(i_a, i_y, i_c1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x123



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x124(so1, io1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: so1, io1
real(kind=8), intent(inout) :: V2(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call g_sigma_ccvv_ccvv_no0_x124(so1, io1, h2_i, x, d2, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_no0_x124



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x124(s_o1, i_o1, V2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_o4, i_o4
! X()  <-- 
! (    1.00000000)  D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_ = X_ &
  + 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * V2_(s_o4, s_o2, s_o3)%array(i_o4, i_o2, i_o3)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x124



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x124(sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x124(sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x124



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x124(s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    2.00000000) X T2(w,y,a,c) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_a,s_c)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * X & 
  * T2_(s_a, s_y, s_w)%array(i_a, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x124



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x125(so1, io1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: so1, io1
real(kind=8), intent(inout) :: V2(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call g_sigma_ccvv_ccvv_no0_x125(so1, io1, h2_i, x, d2, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_no0_x125



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x125(s_o1, i_o1, V2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_o4, i_o4
! X()  <-- 
! (    1.00000000)  D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_ = X_ &
  + 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * V2_(s_o4, s_o2, s_o3)%array(i_o4, i_o2, i_o3)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x125



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x125(sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x125(sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x125



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x125(s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -1.00000000) X T2(y,w,a,c) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_a,s_c)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 1.00000000d+00 & 
  * X & 
  * T2_(s_a, s_w, s_y)%array(i_a, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x125



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x126(so1, io1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: so1, io1
real(kind=8), intent(inout) :: V2(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call g_sigma_ccvv_ccvv_no0_x126(so1, io1, h2_i, x, d2, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_no0_x126



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x126(s_o1, i_o1, V2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_o4, i_o4
! X()  <-- 
! (    1.00000000)  D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_ = X_ &
  + 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * V2_(s_o4, s_o2, s_o3)%array(i_o4, i_o2, i_o3)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x126



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x126(sa, ia, sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x126(sa, ia, sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x126



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x126(s_a, i_a, s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w
! S2(w,y,a,c)  <-- 
! (    2.00000000) X T2(y,w,c,a) 
do s_y = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_c,s_a)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 2.00000000d+00 & 
  * X & 
  * T2_(s_c, s_w, s_y)%array(i_c, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x126



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x127(so1, io1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: so1, io1
real(kind=8), intent(inout) :: V2(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call g_sigma_ccvv_ccvv_no0_x127(so1, io1, h2_i, x, d2, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_no0_x127



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x127(s_o1, i_o1, V2_, X_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_o4, i_o4
! X()  <-- 
! (    1.00000000)  D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_ = X_ &
  + 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * V2_(s_o4, s_o2, s_o3)%array(i_o4, i_o2, i_o3)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x127



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x127(sa, ia, sc, ic, X, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: X
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x127(sa, ia, sc, ic, X, av2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x127



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x127(s_a, i_a, s_c, i_c, X, T2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: X
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y
! S2(w,y,a,c)  <-- 
! (   -1.00000000) X T2(w,y,c,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_c,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 1.00000000d+00 & 
  * X & 
  * T2_(s_c, s_y, s_w)%array(i_c, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x127



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x128(sc, ic, V2, X, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call g_sigma_ccvv_ccvv_no0_x128(sc, ic, h2_i, xv, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xv)

end subroutine g_if_sigma_ccvv_ccvv_no0_x128



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x128(s_c, i_c, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1
! X(c,v1)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c,v1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c,s_v1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
X_(s_v1)%array(i_v1) =  &
    X_(s_v1)%array(i_v1) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_v1)%array(i_o2, i_o1, i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x128



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x128(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x128(sa, ia, sc, ic, av2_i, xv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x128



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x128(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_v1, i_v1
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(y,w,v1,a) X(c,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_v1, s_w, s_y)%array(i_v1, i_w, i_y) & 
  * X_(s_v1)%array(i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x128



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x129(sc, ic, V2, X, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call g_sigma_ccvv_ccvv_no0_x129(sc, ic, h2_i, xv, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xv)

end subroutine g_if_sigma_ccvv_ccvv_no0_x129



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x129(s_c, i_c, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1
! X(c,v1)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c,v1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c,s_v1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
X_(s_v1)%array(i_v1) =  &
    X_(s_v1)%array(i_v1) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_v1)%array(i_o2, i_o1, i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x129



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x129(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x129(sa, ia, sc, ic, av2_i, xv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x129



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x129(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_v1, i_v1
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(w,y,v1,a) X(c,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_v1, s_y, s_w)%array(i_v1, i_y, i_w) & 
  * X_(s_v1)%array(i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x129



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x130(sc, ic, V2, X, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call g_sigma_ccvv_ccvv_no0_x130(sc, ic, h2_i, xv, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xv)

end subroutine g_if_sigma_ccvv_ccvv_no0_x130



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x130(s_c, i_c, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1
! X(c,v1)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c,o1,o2,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c,s_o1) == IEOR(s_o2,s_v1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
X_(s_v1)%array(i_v1) =  &
    X_(s_v1)%array(i_v1) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x130



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x130(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x130(sa, ia, sc, ic, av2_i, xv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x130



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x130(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_v1, i_v1
! S2(w,y,a,c)  <-- 
! (    1.00000000) T2(w,y,v1,a) X(c,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_y, s_w)%array(i_v1, i_y, i_w) & 
  * X_(s_v1)%array(i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x130



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x131(sa, ia, V2, X, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call g_sigma_ccvv_ccvv_no0_x131(sa, ia, h2_i, xv, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xv)

end subroutine g_if_sigma_ccvv_ccvv_no0_x131



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x131(s_a, i_a, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1
! X(a,v1)  <-- 
! (    1.00000000)  D1(o1,o2) V2(a,v1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_v1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
X_(s_v1)%array(i_v1) =  &
    X_(s_v1)%array(i_v1) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_v1)%array(i_o2, i_o1, i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x131



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x131(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x131(sa, ia, sc, ic, av2_i, xv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x131



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x131(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_v1, i_v1
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(w,y,v1,c) X(a,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_v1, s_y, s_w)%array(i_v1, i_y, i_w) & 
  * X_(s_v1)%array(i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x131



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x132(sa, ia, V2, X, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call g_sigma_ccvv_ccvv_no0_x132(sa, ia, h2_i, xv, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xv)

end subroutine g_if_sigma_ccvv_ccvv_no0_x132



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x132(s_a, i_a, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1
! X(a,v1)  <-- 
! (    1.00000000)  D1(o1,o2) V2(a,v1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_v1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
X_(s_v1)%array(i_v1) =  &
    X_(s_v1)%array(i_v1) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_v1)%array(i_o2, i_o1, i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x132



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x132(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x132(sa, ia, sc, ic, av2_i, xv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x132



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x132(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_v1, i_v1
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(y,w,v1,c) X(a,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_v1, s_w, s_y)%array(i_v1, i_w, i_y) & 
  * X_(s_v1)%array(i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x132



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x133(sa, ia, V2, X, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call g_sigma_ccvv_ccvv_no0_x133(sa, ia, h2_i, xv, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xv)

end subroutine g_if_sigma_ccvv_ccvv_no0_x133



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x133(s_a, i_a, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1
! X(a,v1)  <-- 
! (    1.00000000)  D1(o1,o2) V2(a,o1,o2,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_o1) == IEOR(s_o2,s_v1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
X_(s_v1)%array(i_v1) =  &
    X_(s_v1)%array(i_v1) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x133



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x133(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x133(sa, ia, sc, ic, av2_i, xv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x133



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x133(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_v1, i_v1
! S2(w,y,a,c)  <-- 
! (    1.00000000) T2(y,w,v1,c) X(a,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_w, s_y)%array(i_v1, i_w, i_y) & 
  * X_(s_v1)%array(i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x133



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x134(sc, ic, sv1, iv1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: V2(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sc,sv1)

call g_sigma_ccvv_ccvv_no0_x134(sc, ic, sv1, iv1, h2_i, x, d1, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_no0_x134



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x134(s_c, i_c, s_v1, i_v1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X(c,v1)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c,v1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c,s_v1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_v1)%array(i_o2, i_o1, i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x134



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x134(sc, ic, sv1, iv1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), X, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sc,sv1)

call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x134(sc, ic, sv1, iv1, av2_i, x, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x134



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x134(s_c, i_c, s_v1, i_v1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: X_

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(w,y,a,v1) X(c,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) & 
  * X_
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x134



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x135(sc, ic, sv1, iv1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: V2(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sc,sv1)

call g_sigma_ccvv_ccvv_no0_x135(sc, ic, sv1, iv1, h2_i, x, d1, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_no0_x135



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x135(s_c, i_c, s_v1, i_v1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X(c,v1)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c,v1,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c,s_v1) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_o2, s_o1, s_v1)%array(i_o2, i_o1, i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x135



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x135(sc, ic, sv1, iv1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), X, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sc,sv1)

call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x135(sc, ic, sv1, iv1, av2_i, x, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x135



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x135(s_c, i_c, s_v1, i_v1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: X_

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(y,w,a,v1) X(c,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_a, s_w, s_y)%array(i_a, i_w, i_y) & 
  * X_
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x135



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x136(sc, ic, sv1, iv1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: V2(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sc,sv1)

call g_sigma_ccvv_ccvv_no0_x136(sc, ic, sv1, iv1, h2_i, x, d1, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_no0_x136



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x136(s_c, i_c, s_v1, i_v1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X(c,v1)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c,o1,o2,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c,s_o1) == IEOR(s_o2,s_v1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x136



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x136(sc, ic, sv1, iv1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), X, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sc,sv1)

call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x136(sc, ic, sv1, iv1, av2_i, x, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x136



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x136(s_c, i_c, s_v1, i_v1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: X_

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    1.00000000) T2(y,w,a,v1) X(c,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_w, s_y)%array(i_a, i_w, i_y) & 
  * X_
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x136



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x137(sv1, iv1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv1, iv1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv1

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call g_sigma_ccvv_ccvv_no0_x137(sv1, iv1, h2_i, xv, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xv)

end subroutine g_if_sigma_ccvv_ccvv_no0_x137



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x137(s_v1, i_v1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_a, i_a
! X(v1,a)  <-- 
! (    1.00000000)  D1(o1,o2) V2(v1,a,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_v1,s_a) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_v1,s_a) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a)%array(i_a) =  &
    X_(s_a)%array(i_a) &
  + 1.00000000d+00 & 
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

end subroutine g_sigma_ccvv_ccvv_no0_x137



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x137(sc, ic, sv1, iv1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv1

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x137(sc, ic, sv1, iv1, av2_i, xv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x137



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x137(s_c, i_c, s_v1, i_v1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(y,w,c,v1) X(v1,a) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_c,s_v1) .and. &
IEOR(s_v1,s_a) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_c, s_w, s_y)%array(i_c, i_w, i_y) & 
  * X_(s_a)%array(i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x137



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x138(sv1, iv1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv1, iv1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv1

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call g_sigma_ccvv_ccvv_no0_x138(sv1, iv1, h2_i, xv, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xv)

end subroutine g_if_sigma_ccvv_ccvv_no0_x138



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x138(s_v1, i_v1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_a, i_a
! X(v1,a)  <-- 
! (    1.00000000)  D1(o1,o2) V2(v1,a,o1,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_v1,s_a) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_v1,s_a) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a)%array(i_a) =  &
    X_(s_a)%array(i_a) &
  + 1.00000000d+00 & 
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

end subroutine g_sigma_ccvv_ccvv_no0_x138



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x138(sc, ic, sv1, iv1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv1

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x138(sc, ic, sv1, iv1, av2_i, xv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x138



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x138(s_c, i_c, s_v1, i_v1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(w,y,c,v1) X(v1,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_c,s_v1) .and. &
IEOR(s_v1,s_a) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_c, s_y, s_w)%array(i_c, i_y, i_w) & 
  * X_(s_a)%array(i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x138



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x139(sv1, iv1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv1, iv1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv1

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call g_sigma_ccvv_ccvv_no0_x139(sv1, iv1, h2_i, xv, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xv)

end subroutine g_if_sigma_ccvv_ccvv_no0_x139



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x139(s_v1, i_v1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_a, i_a
! X(v1,a)  <-- 
! (    1.00000000)  D1(o1,o2) V2(v1,o2,o1,a) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_v1,s_a) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_v1,s_o2) == IEOR(s_o1,s_a)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a)%array(i_a) =  &
    X_(s_a)%array(i_a) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x139



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x139(sc, ic, sv1, iv1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv1

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x139(sc, ic, sv1, iv1, av2_i, xv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x139



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x139(s_c, i_c, s_v1, i_v1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    1.00000000) T2(w,y,c,v1) X(v1,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_c,s_v1) .and. &
IEOR(s_v1,s_a) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 1.00000000d+00 & 
  * T2_(s_c, s_y, s_w)%array(i_c, i_y, i_w) & 
  * X_(s_a)%array(i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x139



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x140(sc, ic, V2, X, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call g_sigma_ccvv_ccvv_no0_x140(sc, ic, h2_i, xv, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xv)

end subroutine g_if_sigma_ccvv_ccvv_no0_x140



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x140(s_c, i_c, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1
! X(c,v1)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c,o1,o2,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c,s_o1) == IEOR(s_o2,s_v1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
X_(s_v1)%array(i_v1) =  &
    X_(s_v1)%array(i_v1) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x140



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x140(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x140(sa, ia, sc, ic, av2_i, xv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x140



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x140(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_v1, i_v1
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(y,w,v1,a) X(c,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_v1, s_w, s_y)%array(i_v1, i_w, i_y) & 
  * X_(s_v1)%array(i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x140



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x141(sa, ia, V2, X, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call g_sigma_ccvv_ccvv_no0_x141(sa, ia, h2_i, xv, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xv)

end subroutine g_if_sigma_ccvv_ccvv_no0_x141



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x141(s_a, i_a, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1
! X(a,v1)  <-- 
! (    1.00000000)  D1(o1,o2) V2(a,o1,o2,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_a,s_o1) == IEOR(s_o2,s_v1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
X_(s_v1)%array(i_v1) =  &
    X_(s_v1)%array(i_v1) &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x141



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x141(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x141(sa, ia, sc, ic, av2_i, xv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x141



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x141(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_v1, i_v1
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(w,y,v1,c) X(a,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_v1, s_y, s_w)%array(i_v1, i_y, i_w) & 
  * X_(s_v1)%array(i_v1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x141



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x142(sc, ic, sv1, iv1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: V2(*), X
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sc,sv1)

call g_sigma_ccvv_ccvv_no0_x142(sc, ic, sv1, iv1, h2_i, x, d1, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_ccvv_ccvv_no0_x142



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x142(s_c, i_c, s_v1, i_v1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: X_
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2
! X(c,v1)  <-- 
! (    1.00000000)  D1(o1,o2) V2(c,o1,o2,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_c,s_o1) == IEOR(s_o2,s_v1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_ = X_ &
  + 1.00000000d+00 & 
  * D1_(s_o2, s_o1)%array(i_o2, i_o1) & 
  * V2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no0_x142



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x142(sc, ic, sv1, iv1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), X, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sc,sv1)

call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x142(sc, ic, sv1, iv1, av2_i, x, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x142



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x142(s_c, i_c, s_v1, i_v1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: X_

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(w,y,a,v1) X(c,v1) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == 0) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) & 
  * X_
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x142



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x143(sv1, iv1, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv1, iv1
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv1

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call g_sigma_ccvv_ccvv_no0_x143(sv1, iv1, h2_i, xv, d1, nir, nsym, psym)

deallocate(h2_i)
deallocate(xv)

end subroutine g_if_sigma_ccvv_ccvv_no0_x143



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x143(s_v1, i_v1, V2_, X_, D1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_a, i_a
! X(v1,a)  <-- 
! (    1.00000000)  D1(o1,o2) V2(v1,o1,o2,a) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_v1,s_a) == 0 .and. & 
IEOR(s_o1,s_o2) == 0 .and. &
IEOR(s_v1,s_o1) == IEOR(s_o2,s_a)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a)%array(i_a) =  &
    X_(s_a)%array(i_a) &
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

end subroutine g_sigma_ccvv_ccvv_no0_x143



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no1_x143(sc, ic, sv1, iv1, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv1

call set_symblock_xv(sleft, x, nir, nsym, psym) ! -> xv (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no1_x143(sc, ic, sv1, iv1, av2_i, xv, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xv)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no1_x143



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no1_x143(s_c, i_c, s_v1, i_v1, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(y,w,c,v1) X(v1,a) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_c,s_v1) .and. &
IEOR(s_v1,s_a) == 0) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_c, s_w, s_y)%array(i_c, i_w, i_y) & 
  * X_(s_a)%array(i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_ccvv_no1_x143



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x144(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x144(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x144



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x144(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_v2, i_v2, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(w,y,v2,v1) V2(c,v1,a,v2) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_v2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_v2,s_v1) .and. &
IEOR(s_c,s_v1) == IEOR(s_a,s_v2)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v2 = psym(I_BEGIN, I_V, s_v2), psym(I_END, I_V, s_v2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_v2, s_y, s_w)%array(i_v2, i_y, i_w) & 
  * V2_(s_v2, s_a, s_v1)%array(i_v2, i_a, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x144



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x145(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x145(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x145



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x145(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_v2, i_v2, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(y,w,v2,v1) V2(c,v1,a,v2) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_v2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_v2,s_v1) .and. &
IEOR(s_c,s_v1) == IEOR(s_a,s_v2)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v2 = psym(I_BEGIN, I_V, s_v2), psym(I_END, I_V, s_v2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_v2, s_w, s_y)%array(i_v2, i_w, i_y) & 
  * V2_(s_v2, s_a, s_v1)%array(i_v2, i_a, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x145



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x146(sc, ic, sv1, iv1, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x146(sc, ic, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x146



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x146(s_c, i_c, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_v2, i_v2, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    4.00000000) T2(y,w,v2,v1) V2(c,v2,a,v1) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_v2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_y,s_w) == IEOR(s_v2,s_v1) .and. &
IEOR(s_c,s_v2) == IEOR(s_a,s_v1)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v2 = psym(I_BEGIN, I_V, s_v2), psym(I_END, I_V, s_v2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T2_(s_v2, s_w, s_y)%array(i_v2, i_w, i_y) & 
  * V2_(s_v1, s_a, s_v2)%array(i_v1, i_a, i_v2)
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

end subroutine g_sigma_ccvv_ccvv_no0_x146



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_ccvv_no0_x147(sc, ic, sv2, iv2, T2, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv2, iv2
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv2, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_ccvv_no0_x147(sc, ic, sv2, iv2, av2_i, h2_i, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_ccvv_no0_x147



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_ccvv_no0_x147(s_c, i_c, s_v2, i_v2, T2_, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v2, s_v2
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_v1, i_v1, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T2(w,y,v1,v2) V2(c,v1,a,v2) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_v1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_w,s_y) == IEOR(s_v1,s_v2) .and. &
IEOR(s_c,s_v1) == IEOR(s_a,s_v2)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T2_(s_v1, s_y, s_w)%array(i_v1, i_y, i_w) & 
  * V2_(s_v2, s_a, s_v1)%array(i_v2, i_a, i_v1)
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

end subroutine g_sigma_ccvv_ccvv_no0_x147

