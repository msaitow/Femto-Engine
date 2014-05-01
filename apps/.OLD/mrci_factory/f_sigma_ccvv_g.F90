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
subroutine g_if_sigma_ccvv_g_no0_x0(sc, ic, T0, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T0
real(kind=8), intent(inout) :: V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_g_no0_x0(sc, ic, T0, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_g_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_g_no0_x0(s_c, i_c, T0, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: T0
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_y, i_y, s_w, i_w, s_a, i_a
! S2(w,y,a,c)  <-- 
! (    4.00000000) T0 V2(c,y,w,a) 
do s_y = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_y) == IEOR(s_w,s_a)) then
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  + 4.00000000d+00 & 
  * T0 & 
  * V2_(s_a, s_w, s_y)%array(i_a, i_w, i_y)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_g_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_g_no0_x1(sc, ic, T0, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T0
real(kind=8), intent(inout) :: V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_g_no0_x1(sc, ic, T0, h2_i, av2_i2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ccvv_g_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_g_no0_x1(s_c, i_c, T0, V2_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: T0
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_y, i_y, s_a, i_a
! S2(w,y,a,c)  <-- 
! (   -2.00000000) T0 V2(c,w,y,a) 
do s_w = 0, nir-1
do s_y = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_y) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_w) == IEOR(s_y,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_y = psym(I_BEGIN, I_C, s_y), psym(I_END, I_C, s_y)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) =  &
    S2_(s_a, s_y, s_w)%array(i_a, i_y, i_w) &
  - 2.00000000d+00 & 
  * T0 & 
  * V2_(s_a, s_y, s_w)%array(i_a, i_y, i_w)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_g_no0_x1

