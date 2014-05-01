#include "../f_ct.fh"


!                                                              
!   _______________                                  ______    
!  |          |                 .'. .`. `````|`````.~      ~.  
!  |______    |______         .'   `   `.    |    |          | 
!  |          |             .'           `.  |    |          | 
!  |          |___________.'               `.|     `.______.'  
!                                                              



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_g_no0_x0(sa, ia, T0, h, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T0
real(kind=8), intent(inout) :: h(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_g_no0_x0(sa, ia, T0, h1, av2_i2, d2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_g_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_g_no0_x0(s_a, i_a, T0, h_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: T0
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! S2(i,k,m,a)  <-- 
! (    1.00000000) T0 D2(i,m,k,o1) h(o1,a) 
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
  * T0 & 
  * D2_(s_o1, s_k, s_m, s_i)%array(i_o1, i_k, i_m, i_i) & 
  * h_(s_a, s_o1)%array(i_a, i_o1)
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

end subroutine g_sigma_ooov_g_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_g_y1(sc1, ic1, V2, Y0, nir, nsym, psym)

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
call g_sigma_ooov_g_y1(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_ooov_g_y1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_g_y1(s_c1, i_c1, V2_, Y0_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_a, i_a
! Y0(o1,a)  <-- 
! (    1.00000000)  V2(c1,c1,o1,a) 
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o1,s_a) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o1,s_a)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y0_(s_a, s_o1)%array(i_a, i_o1) =  &
    Y0_(s_a, s_o1)%array(i_a, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_o1, s_c1)%array(i_a, i_o1, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_g_y1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_g_no0_x1(sa, ia, T0, Y0, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T0
real(kind=8), intent(inout) :: Y0(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yav(sleft2, Y0, nir, nsym, psym) ! -> Yav (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_g_no0_x1(sa, ia, T0, Yav, av2_i2, d2, nir, nsym, psym)

deallocate(Yav)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_g_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_g_no0_x1(s_a, i_a, T0, Y0_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: T0
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y0_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! S2(i,k,m,a)  <-- 
! (    2.00000000) T0 D2(i,m,k,o1) Y0(o1,a) 
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
  + 2.00000000d+00 & 
  * T0 & 
  * D2_(s_o1, s_k, s_m, s_i)%array(i_o1, i_k, i_m, i_i) & 
  * Y0_(s_a, s_o1)%array(i_a, i_o1)
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

end subroutine g_sigma_ooov_g_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_g_y2(sc1, ic1, V2, Y1, nir, nsym, psym)

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
call g_sigma_ooov_g_y2(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_ooov_g_y2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_g_y2(s_c1, i_c1, V2_, Y1_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_a, i_a
! Y1(o1,a)  <-- 
! (    1.00000000)  V2(c1,o1,c1,a) 
do s_o1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o1,s_a) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_a)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y1_(s_a, s_o1)%array(i_a, i_o1) =  &
    Y1_(s_a, s_o1)%array(i_a, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_c1, s_o1)%array(i_a, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_g_y2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_g_no0_x2(sa, ia, T0, Y1, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T0
real(kind=8), intent(inout) :: Y1(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yav(sleft2, Y1, nir, nsym, psym) ! -> Yav (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_g_no0_x2(sa, ia, T0, Yav, av2_i2, d2, nir, nsym, psym)

deallocate(Yav)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_g_no0_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_g_no0_x2(s_a, i_a, T0, Y1_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: T0
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y1_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! S2(i,k,m,a)  <-- 
! (   -1.00000000) T0 D2(i,m,k,o1) Y1(o1,a) 
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
  - 1.00000000d+00 & 
  * T0 & 
  * D2_(s_o1, s_k, s_m, s_i)%array(i_o1, i_k, i_m, i_i) & 
  * Y1_(s_a, s_o1)%array(i_a, i_o1)
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

end subroutine g_sigma_ooov_g_no0_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_g_no0_x3(sa, ia, T0, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T0
real(kind=8), intent(inout) :: V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_g_no0_x3(sa, ia, T0, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_g_no0_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_g_no0_x3(s_a, i_a, T0, V2_, S2_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: T0
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o2, i_o2, s_o1, i_o1
integer :: s_o3, i_o3
! S2(i,k,m,a)  <-- 
! (    1.00000000) T0 D3(i,m,k,o2,o1,o3) V2(a,o2,o1,o3) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o2,s_o1),s_o3) .and. &
IEOR(s_a,s_o2) == IEOR(s_o1,s_o3)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * T0 & 
  * D3_(s_o3, s_o1, s_o2, s_k, s_m, s_i)%array(i_o3, i_o1, i_o2, i_k, i_m, i_i) & 
  * V2_(s_o3, s_o1, s_o2)%array(i_o3, i_o1, i_o2)
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

end subroutine g_sigma_ooov_g_no0_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_g_no0_x4(sa, ia, T0, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T0
real(kind=8), intent(inout) :: V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_g_no0_x4(sa, ia, T0, h2_i, av2_i2, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_g_no0_x4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ooov_g_no0_x4(s_a, i_a, T0, V2_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: T0
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2, s_m, i_m
! S2(i,k,m,a)  <-- 
! (    1.00000000) T0 D2(i,o1,k,o2) V2(a,o2,m,o1) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_o2) .and. &
IEOR(s_a,s_o2) == IEOR(s_m,s_o1)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    S2_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * T0 & 
  * D2_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i) & 
  * V2_(s_o1, s_m, s_o2)%array(i_o1, i_m, i_o2)
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

end subroutine g_sigma_ooov_g_no0_x4

