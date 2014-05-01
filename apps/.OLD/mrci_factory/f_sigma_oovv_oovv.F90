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
subroutine g_if_sigma_oovv_oovv_y0(h, Y0, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y0(h1, Y0, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_sigma_oovv_oovv_y0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y0(h_, Y0_, nir, nsym, psym)

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

end subroutine g_sigma_oovv_oovv_y0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x0(sc, ic, Y0, T2, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no0_x0(sc, ic, Y0, av2_i, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x0(s_c, i_c, Y0, T2_, S2_, D2_, nir, nsym, psym)

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
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o1, i_o1, s_a, i_a
! S2(i,k,a,c)  <-- 
! (    2.00000000) Y0 D2(i,o2,k,o1) T2(o2,o1,a,c) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o2,s_o1) == IEOR(s_a,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * Y0 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
  * T2_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2)
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

end subroutine g_sigma_oovv_oovv_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y1(sc1, ic1, V2, Y1, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y1(sc1, ic1, h2_i, Y1, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_oovv_oovv_y1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y1(s_c1, i_c1, V2_, Y1_, nir, nsym, psym)

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

end subroutine g_sigma_oovv_oovv_y1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x1(sc, ic, Y1, T2, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no0_x1(sc, ic, Y1, av2_i, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x1(s_c, i_c, Y1, T2_, S2_, D2_, nir, nsym, psym)

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
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o1, i_o1, s_a, i_a
! S2(i,k,a,c)  <-- 
! (    2.00000000) Y1 D2(i,o2,k,o1) T2(o2,o1,a,c) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o2,s_o1) == IEOR(s_a,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * Y1 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
  * T2_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2)
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

end subroutine g_sigma_oovv_oovv_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y2(sc1, ic1, V2, Y2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y2(sc1, ic1, h2_i, Y2, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_oovv_oovv_y2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y2(s_c1, i_c1, V2_, Y2_, nir, nsym, psym)

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

end subroutine g_sigma_oovv_oovv_y2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x2(sc, ic, Y2, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Y2
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no0_x2(sc, ic, Y2, av2_i, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no0_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x2(s_c, i_c, Y2, T2_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y2
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o1, i_o1, s_a, i_a
! S2(i,k,a,c)  <-- 
! (   -1.00000000) Y2 D2(i,o2,k,o1) T2(o2,o1,a,c) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o2,s_o1) == IEOR(s_a,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * Y2 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
  * T2_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2)
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

end subroutine g_sigma_oovv_oovv_no0_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y3(h, Y3, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y3(h1, Y3, nir, nsym, psym)

deallocate(h1)

end subroutine g_if_sigma_oovv_oovv_y3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y3(h_, Y3_, nir, nsym, psym)

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

end subroutine g_sigma_oovv_oovv_y3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x3(sa, ia, sc, ic, Y3, T2, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no0_x3(sa, ia, sc, ic, Y3, av2_i, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no0_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x3(s_a, i_a, s_c, i_c, Y3, T2_, S2_, D2_, nir, nsym, psym)

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
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o1, i_o1
! S2(i,k,a,c)  <-- 
! (    2.00000000) Y3 D2(i,o2,k,o1) T2(o1,o2,c,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_c,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * Y3 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
  * T2_(s_c, s_o2, s_o1)%array(i_c, i_o2, i_o1)
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

end subroutine g_sigma_oovv_oovv_no0_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y4(sc1, ic1, V2, Y4, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y4(sc1, ic1, h2_i, Y4, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_oovv_oovv_y4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y4(s_c1, i_c1, V2_, Y4_, nir, nsym, psym)

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

end subroutine g_sigma_oovv_oovv_y4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x4(sa, ia, sc, ic, Y4, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: Y4
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no0_x4(sa, ia, sc, ic, Y4, av2_i, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no0_x4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x4(s_a, i_a, s_c, i_c, Y4, T2_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y4
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o1, i_o1
! S2(i,k,a,c)  <-- 
! (    2.00000000) Y4 D2(i,o2,k,o1) T2(o1,o2,c,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_c,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * Y4 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
  * T2_(s_c, s_o2, s_o1)%array(i_c, i_o2, i_o1)
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

end subroutine g_sigma_oovv_oovv_no0_x4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y5(sc1, ic1, V2, Y5, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y5(sc1, ic1, h2_i, Y5, nir, nsym, psym)

deallocate(h2_i)

end subroutine g_if_sigma_oovv_oovv_y5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y5(s_c1, i_c1, V2_, Y5_, nir, nsym, psym)

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

end subroutine g_sigma_oovv_oovv_y5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x5(sa, ia, sc, ic, Y5, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: Y5
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no0_x5(sa, ia, sc, ic, Y5, av2_i, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no0_x5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x5(s_a, i_a, s_c, i_c, Y5, T2_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Y5
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o1, i_o1
! S2(i,k,a,c)  <-- 
! (   -1.00000000) Y5 D2(i,o2,k,o1) T2(o1,o2,c,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_c,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * Y5 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
  * T2_(s_c, s_o2, s_o1)%array(i_c, i_o2, i_o1)
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

end subroutine g_sigma_oovv_oovv_no0_x5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x6(h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x6(h1, xaaaa, d3, nir, nsym, psym)

deallocate(h1)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x6(h_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! X(i,o3,k,o2)  <-- 
! (    1.00000000)  D3(i,o3,k,o2,o4,o1) h(o4,o1) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(s_o4,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) =  &
    X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * h_(s_o1, s_o4)%array(i_o1, i_o4)
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

end subroutine g_sigma_oovv_oovv_no0_x6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x6(sc, ic, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x6(sc, ic, av2_i, xaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x6(s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_a, i_a, s_i, i_i, s_k, i_k
! S2(i,k,a,c)  <-- 
! (    1.00000000) T2(o3,o2,a,c) X(i,o3,k,o2) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_a,s_c) .and. &
IEOR(s_i,s_o3) == IEOR(s_k,s_o2)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o2, s_o3)%array(i_a, i_o2, i_o3) & 
  * X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i)
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

end subroutine g_sigma_oovv_oovv_no1_x6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x7(h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x7(h1, xaaaa, d3, nir, nsym, psym)

deallocate(h1)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x7(h_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o3, i_o3, s_o1, i_o1
integer :: s_o4, i_o4
! X(i,o2,k,o3)  <-- 
! (    1.00000000)  D3(i,o2,k,o3,o1,o4) h(o1,o4) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_i,s_o2) == IEOR(s_k,s_o3) .and. & 
IEOR(IEOR(s_i,s_o2),s_k) == IEOR(IEOR(s_o3,s_o1),s_o4) .and. &
IEOR(s_o1,s_o4) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o3, s_k, s_o2, s_i)%array(i_o3, i_k, i_o2, i_i) =  &
    X_(s_o3, s_k, s_o2, s_i)%array(i_o3, i_k, i_o2, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o4, s_o1, s_o3, s_k, s_o2, s_i)%array(i_o4, i_o1, i_o3, i_k, i_o2, i_i) & 
  * h_(s_o4, s_o1)%array(i_o4, i_o1)
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

end subroutine g_sigma_oovv_oovv_no0_x7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x7(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x7(sa, ia, sc, ic, av2_i, xaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x7(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_i, i_i, s_k, i_k
! S2(i,k,a,c)  <-- 
! (    1.00000000) T2(o3,o2,c,a) X(i,o2,k,o3) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_c,s_a) .and. &
IEOR(s_i,s_o2) == IEOR(s_k,s_o3)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * T2_(s_c, s_o2, s_o3)%array(i_c, i_o2, i_o3) & 
  * X_(s_o3, s_k, s_o2, s_i)%array(i_o3, i_k, i_o2, i_i)
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

end subroutine g_sigma_oovv_oovv_no1_x7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x8(sa, ia, sc, ic, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = IEOR(sa,sc)

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call g_sigma_oovv_oovv_no0_x8(sa, ia, sc, ic, av2_i, h1, xaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaa)

end subroutine g_if_sigma_oovv_oovv_no0_x8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x8(s_a, i_a, s_c, i_c, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1
! X(o1,o2,a,c)  <-- 
! (    1.00000000)  T2(o1,o2,v1,a) h(c,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
X_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    X_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1) & 
  * h_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no0_x8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x8(sa, ia, sc, ic, X, S2, nir, nsym, psym)

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

sleft = IEOR(sa,sc)

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x8(sa, ia, sc, ic, xaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x8(s_a, i_a, s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o1, i_o1
! S2(i,k,a,c)  <-- 
! (    1.00000000) D2(i,o2,k,o1) X(o1,o2,a,c) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_a,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
  * X_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_sigma_oovv_oovv_no1_x8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x9(sc, ic, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x9(sc, ic, av2_i, h1, xaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x9(s_c, i_c, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_v1, i_v1, s_a, i_a
! X(o2,o1,c,a)  <-- 
! (    1.00000000)  T2(o2,o1,v1,c) h(a,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_v1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_c,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) =  &
    X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o1, s_o2)%array(i_v1, i_o1, i_o2) & 
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

end subroutine g_sigma_oovv_oovv_no0_x9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x9(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x9(sc, ic, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x9(s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

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
! (    1.00000000) D2(i,o2,k,o1) X(o2,o1,c,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o2,s_o1) == IEOR(s_c,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
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

end subroutine g_sigma_oovv_oovv_no1_x9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x10(sc, ic, sv1, iv1, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x10(sc, ic, sv1, iv1, av2_i, h1, xaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x10(s_c, i_c, s_v1, i_v1, T2_, h_, X_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_a, i_a
! X(o2,o1,a,c)  <-- 
! (    1.00000000)  T2(o2,o1,a,v1) h(c,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_a,s_c) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) =  &
    X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) & 
  * h_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no0_x10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x10(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x10(sc, ic, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x10(s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

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
! (    1.00000000) D2(i,o2,k,o1) X(o2,o1,a,c) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o2,s_o1) == IEOR(s_a,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
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

end subroutine g_sigma_oovv_oovv_no1_x10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x11(sc, ic, sv1, iv1, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x11(sc, ic, sv1, iv1, av2_i, h1, xaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x11(s_c, i_c, s_v1, i_v1, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_a, i_a
! X(o2,o1,c,a)  <-- 
! (    1.00000000)  T2(o2,o1,c,v1) h(a,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_c,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_c,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) =  &
    X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_c, s_o1, s_o2)%array(i_c, i_o1, i_o2) & 
  * h_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no0_x11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x11(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x11(sc, ic, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x11(s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

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

end subroutine g_sigma_oovv_oovv_no1_x11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y12(sc1, ic1, V2, Y6, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y12(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_oovv_oovv_y12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y12(s_c1, i_c1, V2_, Y6_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o4, i_o4
! Y6(o1,o4)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o4) 
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o1,s_o4) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o4)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Y6_(s_o4, s_o1)%array(i_o4, i_o1) =  &
    Y6_(s_o4, s_o1)%array(i_o4, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o4, s_o1, s_c1)%array(i_o4, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_y12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x12(Y6, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: Y6(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y6, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x12(Yaa, xaaaa, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x12(Y6_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y6_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! X(i,o3,k,o2)  <-- 
! (    1.00000000)  D3(i,o3,k,o2,o4,o1) Y6(o1,o4) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(s_o1,s_o4) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) =  &
    X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * Y6_(s_o4, s_o1)%array(i_o4, i_o1)
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

end subroutine g_sigma_oovv_oovv_no0_x12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x12(sc, ic, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x12(sc, ic, av2_i, xaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x12(s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_a, i_a, s_i, i_i, s_k, i_k
! S2(i,k,a,c)  <-- 
! (    2.00000000) T2(o3,o2,a,c) X(i,o3,k,o2) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_a,s_c) .and. &
IEOR(s_i,s_o3) == IEOR(s_k,s_o2)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * T2_(s_a, s_o2, s_o3)%array(i_a, i_o2, i_o3) & 
  * X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i)
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

end subroutine g_sigma_oovv_oovv_no1_x12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y13(sc1, ic1, V2, Y7, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y13(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_oovv_oovv_y13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y13(s_c1, i_c1, V2_, Y7_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o4, i_o4
! Y7(o1,o4)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o4) 
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o1,s_o4) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o4)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Y7_(s_o4, s_o1)%array(i_o4, i_o1) =  &
    Y7_(s_o4, s_o1)%array(i_o4, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o4, s_c1, s_o1)%array(i_o4, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_y13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x13(Y7, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: Y7(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y7, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x13(Yaa, xaaaa, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x13(Y7_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y7_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! X(i,o3,k,o2)  <-- 
! (    1.00000000)  D3(i,o3,k,o2,o4,o1) Y7(o1,o4) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(s_o1,s_o4) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) =  &
    X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * Y7_(s_o4, s_o1)%array(i_o4, i_o1)
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

end subroutine g_sigma_oovv_oovv_no0_x13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x13(sc, ic, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x13(sc, ic, av2_i, xaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x13(s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_a, i_a, s_i, i_i, s_k, i_k
! S2(i,k,a,c)  <-- 
! (   -1.00000000) T2(o3,o2,a,c) X(i,o3,k,o2) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_a,s_c) .and. &
IEOR(s_i,s_o3) == IEOR(s_k,s_o2)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * T2_(s_a, s_o2, s_o3)%array(i_a, i_o2, i_o3) & 
  * X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i)
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

end subroutine g_sigma_oovv_oovv_no1_x13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y14(sc1, ic1, V2, Y8, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y14(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_oovv_oovv_y14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y14(s_c1, i_c1, V2_, Y8_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o4, i_o4
! Y8(o1,o4)  <-- 
! (    1.00000000)  V2(c1,c1,o1,o4) 
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o1,s_o4) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_o4)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Y8_(s_o4, s_o1)%array(i_o4, i_o1) =  &
    Y8_(s_o4, s_o1)%array(i_o4, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o4, s_o1, s_c1)%array(i_o4, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_y14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x14(Y8, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: Y8(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y8, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x14(Yaa, xaaaa, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x14(Y8_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y8_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! X(i,o3,k,o2)  <-- 
! (    1.00000000)  D3(i,o3,k,o2,o4,o1) Y8(o1,o4) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(s_o1,s_o4) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) =  &
    X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * Y8_(s_o4, s_o1)%array(i_o4, i_o1)
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

end subroutine g_sigma_oovv_oovv_no0_x14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x14(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x14(sa, ia, sc, ic, av2_i, xaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x14(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o3, i_o3, s_i, i_i, s_k, i_k
! S2(i,k,a,c)  <-- 
! (    2.00000000) T2(o2,o3,c,a) X(i,o3,k,o2) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_c,s_a) .and. &
IEOR(s_i,s_o3) == IEOR(s_k,s_o2)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * T2_(s_c, s_o3, s_o2)%array(i_c, i_o3, i_o2) & 
  * X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i)
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

end subroutine g_sigma_oovv_oovv_no1_x14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y15(sc1, ic1, V2, Y9, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y15(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_sigma_oovv_oovv_y15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y15(s_c1, i_c1, V2_, Y9_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_o4, i_o4
! Y9(o1,o4)  <-- 
! (    1.00000000)  V2(c1,o1,c1,o4) 
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o1,s_o4) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_o4)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Y9_(s_o4, s_o1)%array(i_o4, i_o1) =  &
    Y9_(s_o4, s_o1)%array(i_o4, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_o4, s_c1, s_o1)%array(i_o4, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_y15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x15(Y9, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: Y9(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y9, nir, nsym, psym) ! -> Yaa (allocate) 
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x15(Yaa, xaaaa, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x15(Y9_, X_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y9_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! X(i,o3,k,o2)  <-- 
! (    1.00000000)  D3(i,o3,k,o2,o4,o1) Y9(o1,o4) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_o3) == IEOR(s_k,s_o2) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(s_o1,s_o4) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) =  &
    X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * Y9_(s_o4, s_o1)%array(i_o4, i_o1)
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

end subroutine g_sigma_oovv_oovv_no0_x15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x15(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x15(sa, ia, sc, ic, av2_i, xaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x15(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o3, i_o3, s_i, i_i, s_k, i_k
! S2(i,k,a,c)  <-- 
! (   -1.00000000) T2(o2,o3,c,a) X(i,o3,k,o2) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_c,s_a) .and. &
IEOR(s_i,s_o3) == IEOR(s_k,s_o2)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * T2_(s_c, s_o3, s_o2)%array(i_c, i_o3, i_o2) & 
  * X_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i)
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

end subroutine g_sigma_oovv_oovv_no1_x15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y16(sc1, ic1, V2, Y10, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y16(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_oovv_oovv_y16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y16(s_c1, i_c1, V2_, Y10_, nir, nsym, psym)

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

integer :: s_c, i_c, s_v1, i_v1
! Y10(c,v1)  <-- 
! (    1.00000000)  V2(c1,c1,c,v1) 
do s_c = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c,s_v1)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y10_(s_v1, s_c)%array(i_v1, i_c) =  &
    Y10_(s_v1, s_c)%array(i_v1, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c, s_c1)%array(i_v1, i_c, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_y16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x16(sa, ia, sc, ic, T2, Y10, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), Y10(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y10, nir, nsym, psym) ! -> Yvv (allocate) 
sleft = IEOR(sa,sc)

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call g_sigma_oovv_oovv_no0_x16(sa, ia, sc, ic, av2_i, Yvv, xaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(xaa)

end subroutine g_if_sigma_oovv_oovv_no0_x16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x16(s_a, i_a, s_c, i_c, T2_, Y10_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y10_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1
! X(o1,o2,a,c)  <-- 
! (    1.00000000)  T2(o1,o2,v1,a) Y10(c,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
X_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    X_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1) & 
  * Y10_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no0_x16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x16(sa, ia, sc, ic, X, S2, nir, nsym, psym)

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

sleft = IEOR(sa,sc)

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x16(sa, ia, sc, ic, xaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x16(s_a, i_a, s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o1, i_o1
! S2(i,k,a,c)  <-- 
! (    2.00000000) D2(i,o2,k,o1) X(o1,o2,a,c) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_a,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
  * X_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_sigma_oovv_oovv_no1_x16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y17(sc1, ic1, V2, Y11, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y17(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_oovv_oovv_y17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y17(s_c1, i_c1, V2_, Y11_, nir, nsym, psym)

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

integer :: s_c, i_c, s_v1, i_v1
! Y11(c,v1)  <-- 
! (    1.00000000)  V2(c1,c,c1,v1) 
do s_c = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_c1,s_c) == IEOR(s_c1,s_v1)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y11_(s_v1, s_c)%array(i_v1, i_c) =  &
    Y11_(s_v1, s_c)%array(i_v1, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_c)%array(i_v1, i_c1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_y17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x17(sa, ia, sc, ic, T2, Y11, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), Y11(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y11, nir, nsym, psym) ! -> Yvv (allocate) 
sleft = IEOR(sa,sc)

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call g_sigma_oovv_oovv_no0_x17(sa, ia, sc, ic, av2_i, Yvv, xaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(xaa)

end subroutine g_if_sigma_oovv_oovv_no0_x17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x17(s_a, i_a, s_c, i_c, T2_, Y11_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y11_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1
! X(o1,o2,a,c)  <-- 
! (    1.00000000)  T2(o1,o2,v1,a) Y11(c,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
X_(s_o2, s_o1)%array(i_o2, i_o1) =  &
    X_(s_o2, s_o1)%array(i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1) & 
  * Y11_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no0_x17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x17(sa, ia, sc, ic, X, S2, nir, nsym, psym)

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

sleft = IEOR(sa,sc)

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x17(sa, ia, sc, ic, xaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x17(s_a, i_a, s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o1, i_o1
! S2(i,k,a,c)  <-- 
! (   -1.00000000) D2(i,o2,k,o1) X(o1,o2,a,c) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_a,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
  * X_(s_o2, s_o1)%array(i_o2, i_o1)
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

end subroutine g_sigma_oovv_oovv_no1_x17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y18(sc1, ic1, V2, Y12, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y18(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_oovv_oovv_y18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y18(s_c1, i_c1, V2_, Y12_, nir, nsym, psym)

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

end subroutine g_sigma_oovv_oovv_y18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x18(sc, ic, T2, Y12, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), Y12(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y12, nir, nsym, psym) ! -> Yvv (allocate) 
sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x18(sc, ic, av2_i, Yvv, xaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x18(s_c, i_c, T2_, Y12_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y12_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_v1, i_v1, s_a, i_a
! X(o2,o1,c,a)  <-- 
! (    1.00000000)  T2(o2,o1,v1,c) Y12(a,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_v1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_c,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) =  &
    X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o1, s_o2)%array(i_v1, i_o1, i_o2) & 
  * Y12_(s_v1, s_a)%array(i_v1, i_a)
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

end subroutine g_sigma_oovv_oovv_no0_x18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x18(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x18(sc, ic, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x18(s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

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
! (    2.00000000) D2(i,o2,k,o1) X(o2,o1,c,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o2,s_o1) == IEOR(s_c,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
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

end subroutine g_sigma_oovv_oovv_no1_x18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y19(sc1, ic1, V2, Y13, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y19(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_oovv_oovv_y19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y19(s_c1, i_c1, V2_, Y13_, nir, nsym, psym)

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

end subroutine g_sigma_oovv_oovv_y19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x19(sc, ic, T2, Y13, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), Y13(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y13, nir, nsym, psym) ! -> Yvv (allocate) 
sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x19(sc, ic, av2_i, Yvv, xaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x19(s_c, i_c, T2_, Y13_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y13_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_v1, i_v1, s_a, i_a
! X(o2,o1,c,a)  <-- 
! (    1.00000000)  T2(o2,o1,v1,c) Y13(a,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_v1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_c,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) =  &
    X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o1, s_o2)%array(i_v1, i_o1, i_o2) & 
  * Y13_(s_v1, s_a)%array(i_v1, i_a)
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

end subroutine g_sigma_oovv_oovv_no0_x19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x19(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x19(sc, ic, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x19(s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

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
! (   -1.00000000) D2(i,o2,k,o1) X(o2,o1,c,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o2,s_o1) == IEOR(s_c,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
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

end subroutine g_sigma_oovv_oovv_no1_x19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y20(sc1, ic1, V2, Y14, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y20(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_oovv_oovv_y20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y20(s_c1, i_c1, V2_, Y14_, nir, nsym, psym)

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

integer :: s_c, i_c, s_v1, i_v1
! Y14(c,v1)  <-- 
! (    1.00000000)  V2(c1,c1,c,v1) 
do s_c = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_c,s_v1)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y14_(s_v1, s_c)%array(i_v1, i_c) =  &
    Y14_(s_v1, s_c)%array(i_v1, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c, s_c1)%array(i_v1, i_c, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_y20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x20(sc, ic, sv1, iv1, T2, Y14, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y14(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y14, nir, nsym, psym) ! -> Yvv (allocate) 
sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x20(sc, ic, sv1, iv1, av2_i, Yvv, xaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x20(s_c, i_c, s_v1, i_v1, T2_, Y14_, X_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: Y14_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_a, i_a
! X(o2,o1,a,c)  <-- 
! (    1.00000000)  T2(o2,o1,a,v1) Y14(c,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_a,s_c) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) =  &
    X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) & 
  * Y14_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no0_x20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x20(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x20(sc, ic, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x20(s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

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
! (    2.00000000) D2(i,o2,k,o1) X(o2,o1,a,c) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o2,s_o1) == IEOR(s_a,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
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

end subroutine g_sigma_oovv_oovv_no1_x20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y21(sc1, ic1, V2, Y15, nir, nsym, psym)

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
call g_sigma_oovv_oovv_y21(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_oovv_oovv_y21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y21(s_c1, i_c1, V2_, Y15_, nir, nsym, psym)

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

integer :: s_c, i_c, s_v1, i_v1
! Y15(c,v1)  <-- 
! (    1.00000000)  V2(c1,c,c1,v1) 
do s_c = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_c,s_v1) == 0 .and. & 
IEOR(s_c1,s_c) == IEOR(s_c1,s_v1)) then
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y15_(s_v1, s_c)%array(i_v1, i_c) =  &
    Y15_(s_v1, s_c)%array(i_v1, i_c) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_c)%array(i_v1, i_c1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_y21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x21(sc, ic, sv1, iv1, T2, Y15, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y15(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y15, nir, nsym, psym) ! -> Yvv (allocate) 
sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x21(sc, ic, sv1, iv1, av2_i, Yvv, xaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x21(s_c, i_c, s_v1, i_v1, T2_, Y15_, X_, nir, nsym, psym)

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
type(symblock2), intent(inout) :: Y15_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_a, i_a
! X(o2,o1,a,c)  <-- 
! (    1.00000000)  T2(o2,o1,a,v1) Y15(c,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_a,s_c) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) =  &
    X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) & 
  * Y15_(s_v1, s_c)%array(i_v1, i_c)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no0_x21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x21(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x21(sc, ic, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x21(s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

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
! (   -1.00000000) D2(i,o2,k,o1) X(o2,o1,a,c) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o2,s_o1) == IEOR(s_a,s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
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

end subroutine g_sigma_oovv_oovv_no1_x21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y22(sc1, ic1, V2, Y16, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y16, nir, nsym, psym) ! -> Yvv (allocate) 
call g_sigma_oovv_oovv_y22(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_oovv_oovv_y22



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y22(s_c1, i_c1, V2_, Y16_, nir, nsym, psym)

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

integer :: s_a, i_a, s_v1, i_v1
! Y16(a,v1)  <-- 
! (    1.00000000)  V2(c1,c1,a,v1) 
do s_a = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_v1)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y16_(s_v1, s_a)%array(i_v1, i_a) =  &
    Y16_(s_v1, s_a)%array(i_v1, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_a, s_c1)%array(i_v1, i_a, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_y22



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x22(sc, ic, sv1, iv1, T2, Y16, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y16(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y16, nir, nsym, psym) ! -> Yvv (allocate) 
sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x22(sc, ic, sv1, iv1, av2_i, Yvv, xaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x22



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x22(s_c, i_c, s_v1, i_v1, T2_, Y16_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y16_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_a, i_a
! X(o1,o2,c,a)  <-- 
! (    1.00000000)  T2(o1,o2,c,v1) Y16(a,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_c,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_c,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1) =  &
    X_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_c, s_o2, s_o1)%array(i_c, i_o2, i_o1) & 
  * Y16_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no0_x22



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x22(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x22(sc, ic, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x22



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x22(s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

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
! (    2.00000000) D2(i,o2,k,o1) X(o1,o2,c,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_c,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 2.00000000d+00 & 
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

end subroutine g_sigma_oovv_oovv_no1_x22



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_y23(sc1, ic1, V2, Y17, nir, nsym, psym)

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
call set_symblock_Yvv(sleft2, Y17, nir, nsym, psym) ! -> Yvv (allocate) 
call g_sigma_oovv_oovv_y23(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_sigma_oovv_oovv_y23



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_y23(s_c1, i_c1, V2_, Y17_, nir, nsym, psym)

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

integer :: s_a, i_a, s_v1, i_v1
! Y17(a,v1)  <-- 
! (    1.00000000)  V2(c1,a,c1,v1) 
do s_a = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_a,s_v1) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_v1)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y17_(s_v1, s_a)%array(i_v1, i_a) =  &
    Y17_(s_v1, s_a)%array(i_v1, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_a)%array(i_v1, i_c1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_y23



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x23(sc, ic, sv1, iv1, T2, Y17, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y17(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y17, nir, nsym, psym) ! -> Yvv (allocate) 
sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x23(sc, ic, sv1, iv1, av2_i, Yvv, xaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yvv)
deallocate(xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x23



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x23(s_c, i_c, s_v1, i_v1, T2_, Y17_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y17_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_a, i_a
! X(o1,o2,c,a)  <-- 
! (    1.00000000)  T2(o1,o2,c,v1) Y17(a,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_c,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_c,s_v1) .and. &
IEOR(s_a,s_v1) == 0) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1) =  &
    X_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_c, s_o2, s_o1)%array(i_c, i_o2, i_o1) & 
  * Y17_(s_v1, s_a)%array(i_v1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no0_x23



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x23(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x23(sc, ic, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x23



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x23(s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

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
! (   -1.00000000) D2(i,o2,k,o1) X(o1,o2,c,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o1,s_o2) == IEOR(s_c,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  - 1.00000000d+00 & 
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

end subroutine g_sigma_oovv_oovv_no1_x23



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x24(so1, io1, so5, io5, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: so1, io1, so5, io5
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x24(so1, io1, so5, io5, h2_i, xaaaa, d4_ij, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x24



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x24(s_o1, i_o1, s_o5, i_o5, V2_, X_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
integer, intent(in) :: i_o5, s_o5
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o4, i_o4, s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2
integer :: s_o6, i_o6
! X(o4,i,o3,k)  <-- 
! (    1.00000000)  D4(o1,o5,o4,i,o3,k,o2,o6) V2(o1,o5,o2,o6) 
do s_o4 = 0, nir-1
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o6 = 0, nir-1
if( &
IEOR(s_o4,s_i) == IEOR(s_o3,s_k) .and. & 
IEOR(IEOR(s_o1,s_o5),IEOR(s_o4,s_i)) == IEOR(IEOR(s_o3,s_k),IEOR(s_o2,s_o6)) .and. &
IEOR(s_o1,s_o5) == IEOR(s_o2,s_o6)) then
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o6 = psym(I_BEGIN, I_O, s_o6), psym(I_END, I_O, s_o6)
X_(s_k, s_o3, s_i, s_o4)%array(i_k, i_o3, i_i, i_o4) =  &
    X_(s_k, s_o3, s_i, s_o4)%array(i_k, i_o3, i_i, i_o4) &
  + 1.00000000d+00 & 
  * D4_(s_o6, s_o2, s_k, s_o3, s_i, s_o4)%array(i_o6, i_o2, i_k, i_o3, i_i, i_o4) & 
  * V2_(s_o6, s_o2, s_o5)%array(i_o6, i_o2, i_o5)
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

end subroutine g_sigma_oovv_oovv_no0_x24



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x24(sc, ic, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x24(sc, ic, av2_i, xaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x24



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x24(s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o4, i_o4, s_o3, i_o3, s_a, i_a, s_i, i_i, s_k, i_k
! S2(i,k,a,c)  <-- 
! (    0.50000000) T2(o4,o3,a,c) X(o4,i,o3,k) 
do s_o4 = 0, nir-1
do s_o3 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o4,s_o3) == IEOR(s_a,s_c) .and. &
IEOR(s_o4,s_i) == IEOR(s_o3,s_k)) then
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 0.50000000d+00 & 
  * T2_(s_a, s_o3, s_o4)%array(i_a, i_o3, i_o4) & 
  * X_(s_k, s_o3, s_i, s_o4)%array(i_k, i_o3, i_i, i_o4)
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

end subroutine g_sigma_oovv_oovv_no1_x24



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x25(so1, io1, so5, io5, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: so1, io1, so5, io5
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x25(so1, io1, so5, io5, h2_i, xaaaa, d4_ij, nir, nsym, psym)

deallocate(h2_i)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x25



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x25(s_o1, i_o1, s_o5, i_o5, V2_, X_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_o1, s_o1
integer, intent(in) :: i_o5, s_o5
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o4, i_o4, s_o2, i_o2
integer :: s_o6, i_o6
! X(i,o3,k,o4)  <-- 
! (    1.00000000)  D4(o1,o5,i,o3,k,o4,o2,o6) V2(o1,o5,o2,o6) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
do s_o6 = 0, nir-1
if( &
IEOR(s_i,s_o3) == IEOR(s_k,s_o4) .and. & 
IEOR(IEOR(s_o1,s_o5),IEOR(s_i,s_o3)) == IEOR(IEOR(s_k,s_o4),IEOR(s_o2,s_o6)) .and. &
IEOR(s_o1,s_o5) == IEOR(s_o2,s_o6)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o6 = psym(I_BEGIN, I_O, s_o6), psym(I_END, I_O, s_o6)
X_(s_o4, s_k, s_o3, s_i)%array(i_o4, i_k, i_o3, i_i) =  &
    X_(s_o4, s_k, s_o3, s_i)%array(i_o4, i_k, i_o3, i_i) &
  + 1.00000000d+00 & 
  * D4_(s_o6, s_o2, s_o4, s_k, s_o3, s_i)%array(i_o6, i_o2, i_o4, i_k, i_o3, i_i) & 
  * V2_(s_o6, s_o2, s_o5)%array(i_o6, i_o2, i_o5)
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

end subroutine g_sigma_oovv_oovv_no0_x25



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x25(sa, ia, sc, ic, T2, X, S2, nir, nsym, psym)

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
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x25(sa, ia, sc, ic, av2_i, xaaaa, av2_i2, nir, nsym, psym)

deallocate(av2_i)
deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x25



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x25(s_a, i_a, s_c, i_c, T2_, X_, S2_, nir, nsym, psym)

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
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o4, i_o4, s_o3, i_o3, s_i, i_i, s_k, i_k
! S2(i,k,a,c)  <-- 
! (    0.50000000) T2(o4,o3,c,a) X(i,o3,k,o4) 
do s_o4 = 0, nir-1
do s_o3 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o4,s_o3) == IEOR(s_c,s_a) .and. &
IEOR(s_i,s_o3) == IEOR(s_k,s_o4)) then
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 0.50000000d+00 & 
  * T2_(s_c, s_o3, s_o4)%array(i_c, i_o3, i_o4) & 
  * X_(s_o4, s_k, s_o3, s_i)%array(i_o4, i_k, i_o3, i_i)
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

end subroutine g_sigma_oovv_oovv_no1_x25



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x26(sa, ia, sc, ic, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa,sc)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x26(sa, ia, sc, ic, av2_i, h2_i, xaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x26



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x26(s_a, i_a, s_c, i_c, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o3, i_o3, s_v1, i_v1, s_o1, i_o1, s_o4, i_o4
! X(o2,o3,o1,o4,a,c)  <-- 
! (    1.00000000)  T2(o2,o3,v1,a) V2(c,v1,o1,o4) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_v1 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(IEOR(s_o2,s_o3),s_o1) == IEOR(IEOR(s_o4,s_a),s_c) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_v1) == IEOR(s_o1,s_o4)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o1, s_o3, s_o2)%array(i_o4, i_o1, i_o3, i_o2) =  &
    X_(s_o4, s_o1, s_o3, s_o2)%array(i_o4, i_o1, i_o3, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o3, s_o2)%array(i_v1, i_o3, i_o2) & 
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

end subroutine g_sigma_oovv_oovv_no0_x26



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x26(sa, ia, sc, ic, X, S2, nir, nsym, psym)

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

sleft = IEOR(sa,sc)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x26(sa, ia, sc, ic, xaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x26



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x26(s_a, i_a, s_c, i_c, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! S2(i,k,a,c)  <-- 
! (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o2,o3,o1,o4,a,c) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(IEOR(s_o2,s_o3),s_o1) == IEOR(IEOR(s_o4,s_a),s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * X_(s_o4, s_o1, s_o3, s_o2)%array(i_o4, i_o1, i_o3, i_o2)
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

end subroutine g_sigma_oovv_oovv_no1_x26



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x27(sa, ia, sc, ic, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sc,sa)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x27(sa, ia, sc, ic, av2_i, h2_i, xaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x27



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x27(s_a, i_a, s_c, i_c, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_v1, i_v1, s_o1, i_o1, s_o4, i_o4
! X(o3,o2,o1,o4,c,a)  <-- 
! (    1.00000000)  T2(o3,o2,v1,c) V2(a,v1,o1,o4) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(IEOR(s_o3,s_o2),s_o1) == IEOR(IEOR(s_o4,s_c),s_a) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_v1) == IEOR(s_o1,s_o4)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o1, s_o2, s_o3)%array(i_o4, i_o1, i_o2, i_o3) =  &
    X_(s_o4, s_o1, s_o2, s_o3)%array(i_o4, i_o1, i_o2, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o2, s_o3)%array(i_v1, i_o2, i_o3) & 
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

end subroutine g_sigma_oovv_oovv_no0_x27



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x27(sa, ia, sc, ic, X, S2, nir, nsym, psym)

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

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x27(sa, ia, sc, ic, xaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x27



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x27(s_a, i_a, s_c, i_c, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! S2(i,k,a,c)  <-- 
! (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o3,o2,o1,o4,c,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(IEOR(s_o3,s_o2),s_o1) == IEOR(IEOR(s_o4,s_c),s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * X_(s_o4, s_o1, s_o2, s_o3)%array(i_o4, i_o1, i_o2, i_o3)
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

end subroutine g_sigma_oovv_oovv_no1_x27



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x28(sc, ic, sv1, iv1, T2, V2, X, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no0_x28(sc, ic, sv1, iv1, av2_i, h2_i, xaaaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaav)

end subroutine g_if_sigma_oovv_oovv_no0_x28



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x28(s_c, i_c, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

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

integer :: s_o3, i_o3, s_o2, i_o2, s_a, i_a, s_o1, i_o1, s_o4, i_o4
! X(o3,o2,o1,o4,a,c)  <-- 
! (    1.00000000)  T2(o3,o2,a,v1) V2(c,v1,o1,o4) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(IEOR(s_o3,s_o2),s_o1) == IEOR(IEOR(s_o4,s_a),s_c) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_v1) == IEOR(s_o1,s_o4)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_a, s_o4, s_o1, s_o2, s_o3)%array(i_a, i_o4, i_o1, i_o2, i_o3) =  &
    X_(s_a, s_o4, s_o1, s_o2, s_o3)%array(i_a, i_o4, i_o1, i_o2, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o2, s_o3)%array(i_a, i_o2, i_o3) & 
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

end subroutine g_sigma_oovv_oovv_no0_x28



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x28(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x28(sc, ic, xaaaav, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x28



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x28(s_c, i_c, X_, S2_, D3_, nir, nsym, psym)

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
! (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o3,o2,o1,o4,a,c) 
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
IEOR(IEOR(s_o3,s_o2),s_o1) == IEOR(IEOR(s_o4,s_a),s_c)) then
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
  * X_(s_a, s_o4, s_o1, s_o2, s_o3)%array(i_a, i_o4, i_o1, i_o2, i_o3)
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

end subroutine g_sigma_oovv_oovv_no1_x28



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x29(sc, ic, sv1, iv1, T2, V2, X, nir, nsym, psym)

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
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc

call set_symblock_xaaaav(sleft, x, nir, nsym, psym) ! -> xaaaav (allocate) 
call g_sigma_oovv_oovv_no0_x29(sc, ic, sv1, iv1, av2_i, h2_i, xaaaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaav)

end subroutine g_if_sigma_oovv_oovv_no0_x29



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x29(s_c, i_c, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o3, i_o3, s_a, i_a, s_o1, i_o1, s_o4, i_o4
! X(o2,o3,o1,o4,c,a)  <-- 
! (    1.00000000)  T2(o2,o3,c,v1) V2(v1,a,o1,o4) 
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_a = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(IEOR(s_o2,s_o3),s_o1) == IEOR(IEOR(s_o4,s_c),s_a) .and. & 
IEOR(s_o2,s_o3) == IEOR(s_c,s_v1) .and. &
IEOR(s_v1,s_a) == IEOR(s_o1,s_o4)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_a, s_o4, s_o1, s_o3, s_o2)%array(i_a, i_o4, i_o1, i_o3, i_o2) =  &
    X_(s_a, s_o4, s_o1, s_o3, s_o2)%array(i_a, i_o4, i_o1, i_o3, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_c, s_o3, s_o2)%array(i_c, i_o3, i_o2) & 
  * V2_(s_o4, s_o1, s_a)%array(i_o4, i_o1, i_a)
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

end subroutine g_sigma_oovv_oovv_no0_x29



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x29(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x29(sc, ic, xaaaav, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x29



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x29(s_c, i_c, X_, S2_, D3_, nir, nsym, psym)

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
! (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o2,o3,o1,o4,c,a) 
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
IEOR(IEOR(s_o2,s_o3),s_o1) == IEOR(IEOR(s_o4,s_c),s_a)) then
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
  * X_(s_a, s_o4, s_o1, s_o3, s_o2)%array(i_a, i_o4, i_o1, i_o3, i_o2)
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

end subroutine g_sigma_oovv_oovv_no1_x29



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x30(sa, ia, sc, ic, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa,sc)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x30(sa, ia, sc, ic, av2_i, h2_i, xaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x30



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x30(s_a, i_a, s_c, i_c, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3, s_v1, i_v1, s_o2, i_o2, s_o4, i_o4
! X(o1,o3,o2,o4,a,c)  <-- 
! (    1.00000000)  T2(o1,o3,v1,a) V2(c,o2,o4,v1) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_v1 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(IEOR(s_o1,s_o3),s_o2) == IEOR(IEOR(s_o4,s_a),s_c) .and. & 
IEOR(s_o1,s_o3) == IEOR(s_v1,s_a) .and. &
IEOR(s_c,s_o2) == IEOR(s_o4,s_v1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) =  &
    X_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o3, s_o1)%array(i_v1, i_o3, i_o1) & 
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

end subroutine g_sigma_oovv_oovv_no0_x30



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x30(sa, ia, sc, ic, X, S2, nir, nsym, psym)

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

sleft = IEOR(sa,sc)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x30(sa, ia, sc, ic, xaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x30



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x30(s_a, i_a, s_c, i_c, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! S2(i,k,a,c)  <-- 
! (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o1,o3,o2,o4,a,c) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(IEOR(s_o1,s_o3),s_o2) == IEOR(IEOR(s_o4,s_a),s_c)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * X_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1)
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

end subroutine g_sigma_oovv_oovv_no1_x30



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x31(sa, ia, sc, ic, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sc,sa)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x31(sa, ia, sc, ic, av2_i, h2_i, xaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x31



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x31(s_a, i_a, s_c, i_c, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1, s_o3, i_o3, s_o4, i_o4
! X(o1,o2,o3,o4,c,a)  <-- 
! (    1.00000000)  T2(o1,o2,v1,c) V2(a,o3,o4,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(IEOR(s_o1,s_o2),s_o3) == IEOR(IEOR(s_o4,s_c),s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_v1,s_c) .and. &
IEOR(s_a,s_o3) == IEOR(s_o4,s_v1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_o4, s_o3, s_o2, s_o1)%array(i_o4, i_o3, i_o2, i_o1) =  &
    X_(s_o4, s_o3, s_o2, s_o1)%array(i_o4, i_o3, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1) & 
  * V2_(s_v1, s_o4, s_o3)%array(i_v1, i_o4, i_o3)
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

end subroutine g_sigma_oovv_oovv_no0_x31



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x31(sa, ia, sc, ic, X, S2, nir, nsym, psym)

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

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x31(sa, ia, sc, ic, xaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x31



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x31(s_a, i_a, s_c, i_c, X_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! S2(i,k,a,c)  <-- 
! (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o1,o2,o3,o4,c,a) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(IEOR(s_o1,s_o2),s_o3) == IEOR(IEOR(s_o4,s_c),s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i) & 
  * X_(s_o4, s_o3, s_o2, s_o1)%array(i_o4, i_o3, i_o2, i_o1)
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

end subroutine g_sigma_oovv_oovv_no1_x31



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x32(sc, ic, sv1, iv1, T2, V2, X, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no0_x32(sc, ic, sv1, iv1, av2_i, h2_i, xaaaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaav)

end subroutine g_if_sigma_oovv_oovv_no0_x32



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x32(s_c, i_c, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

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

integer :: s_o3, i_o3, s_o1, i_o1, s_a, i_a, s_o2, i_o2, s_o4, i_o4
! X(o3,o1,o2,o4,a,c)  <-- 
! (    1.00000000)  T2(o3,o1,a,v1) V2(c,o2,o4,v1) 
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(IEOR(s_o3,s_o1),s_o2) == IEOR(IEOR(s_o4,s_a),s_c) .and. & 
IEOR(s_o3,s_o1) == IEOR(s_a,s_v1) .and. &
IEOR(s_c,s_o2) == IEOR(s_o4,s_v1)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
X_(s_a, s_o4, s_o2, s_o1, s_o3)%array(i_a, i_o4, i_o2, i_o1, i_o3) =  &
    X_(s_a, s_o4, s_o2, s_o1, s_o3)%array(i_a, i_o4, i_o2, i_o1, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_a, s_o1, s_o3)%array(i_a, i_o1, i_o3) & 
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

end subroutine g_sigma_oovv_oovv_no0_x32



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x32(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x32(sc, ic, xaaaav, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x32



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x32(s_c, i_c, X_, S2_, D3_, nir, nsym, psym)

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
! (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o3,o1,o2,o4,a,c) 
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
IEOR(IEOR(s_o3,s_o1),s_o2) == IEOR(IEOR(s_o4,s_a),s_c)) then
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
  * X_(s_a, s_o4, s_o2, s_o1, s_o3)%array(i_a, i_o4, i_o2, i_o1, i_o3)
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

end subroutine g_sigma_oovv_oovv_no1_x32



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x33(sc, ic, sv1, iv1, T2, V2, X, nir, nsym, psym)

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
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc

call set_symblock_xaaaav(sleft, x, nir, nsym, psym) ! -> xaaaav (allocate) 
call g_sigma_oovv_oovv_no0_x33(sc, ic, sv1, iv1, av2_i, h2_i, xaaaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaav)

end subroutine g_if_sigma_oovv_oovv_no0_x33



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x33(s_c, i_c, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o4, i_o4, s_o1, i_o1, s_o2, i_o2, s_a, i_a
! X(o3,o4,o1,o2,c,a)  <-- 
! (    1.00000000)  T2(o3,o4,c,v1) V2(v1,o1,o2,a) 
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(IEOR(s_o3,s_o4),s_o1) == IEOR(IEOR(s_o2,s_c),s_a) .and. & 
IEOR(s_o3,s_o4) == IEOR(s_c,s_v1) .and. &
IEOR(s_v1,s_o1) == IEOR(s_o2,s_a)) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o2, s_o1, s_o4, s_o3)%array(i_a, i_o2, i_o1, i_o4, i_o3) =  &
    X_(s_a, s_o2, s_o1, s_o4, s_o3)%array(i_a, i_o2, i_o1, i_o4, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_c, s_o4, s_o3)%array(i_c, i_o4, i_o3) & 
  * V2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
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

end subroutine g_sigma_oovv_oovv_no0_x33



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x33(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x33(sc, ic, xaaaav, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x33



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x33(s_c, i_c, X_, S2_, D3_, nir, nsym, psym)

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
! (    1.00000000) D3(i,o2,k,o3,o1,o4) X(o3,o4,o1,o2,c,a) 
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
IEOR(IEOR(s_o3,s_o4),s_o1) == IEOR(IEOR(s_o2,s_c),s_a)) then
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
  * X_(s_a, s_o2, s_o1, s_o4, s_o3)%array(i_a, i_o2, i_o1, i_o4, i_o3)
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

end subroutine g_sigma_oovv_oovv_no1_x33



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x34(sc, ic, sv1, iv1, T2, V2, X, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no0_x34(sc, ic, sv1, iv1, av2_i, h2_i, xaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x34



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x34(s_c, i_c, s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

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

integer :: s_o2, i_o2, s_o1, i_o1, s_v2, i_v2, s_a, i_a
! X(o2,o1,c,a)  <-- 
! (    1.00000000)  T2(o2,o1,v2,v1) V2(c,v1,a,v2) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_v2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_c,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_v2,s_v1) .and. &
IEOR(s_c,s_v1) == IEOR(s_a,s_v2)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_v2 = psym(I_BEGIN, I_V, s_v2), psym(I_END, I_V, s_v2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) =  &
    X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_v2, s_o1, s_o2)%array(i_v2, i_o1, i_o2) & 
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

end subroutine g_sigma_oovv_oovv_no0_x34



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x34(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x34(sc, ic, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x34



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x34(s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

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
! (    1.00000000) D2(i,o2,k,o1) X(o2,o1,c,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_o2,s_o1) == IEOR(s_c,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
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

end subroutine g_sigma_oovv_oovv_no1_x34



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x35(sc, ic, sv2, iv2, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv2, iv2
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv2, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x35(sc, ic, sv2, iv2, av2_i, h2_i, xaav, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x35



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x35(s_c, i_c, s_v2, i_v2, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v2, s_v2
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_v1, i_v1, s_a, i_a
! X(o2,o1,c,a)  <-- 
! (    1.00000000)  T2(o2,o1,v1,v2) V2(c,v1,a,v2) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_v1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_c,s_a) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_v1,s_v2) .and. &
IEOR(s_c,s_v1) == IEOR(s_a,s_v2)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) =  &
    X_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o1, s_o2)%array(i_v1, i_o1, i_o2) & 
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

end subroutine g_sigma_oovv_oovv_no0_x35



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x35(sc, ic, X, S2, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no1_x35(sc, ic, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x35



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x35(s_c, i_c, X_, S2_, D2_, nir, nsym, psym)

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

end subroutine g_sigma_oovv_oovv_no1_x35

