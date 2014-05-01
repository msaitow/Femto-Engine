#include <sci/icmr/fsrc/f_mr.fh>


!  8888888888                     888                  
!  888                            888                  
!  888                            888                  
!  8888888  .d88b.  88888b.d88b.  888888  .d88b.       
!  888     d8P  Y8b 888 "888 "88b 888    d88""88b  
!  888     88888888 888  888  888 888    888  888      
!  888     Y8b.     888  888  888 Y88b.  Y88..88P      
!  888      "Y8888  888  888  888  "Y888  "Y88P"   

!                                    Generated date : Sun Apr 20 10:26:26 2014



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_g_no0_x0_type1_eri_v &
  (sb, ib, T0, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: T0
real(kind=8), intent(inout) :: V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_g_no0_x0_type1_eri_v &
  (sb, ib, T0, h2_i, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)

end subroutine g_if_sigma_ccvv_g_no0_x0_type1_eri_v



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_g_no0_x0_type1_eri_v &
  (s_b, i_b, T0, V2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: T0
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a, i_a
! S2(w,x,a,b) += (    4.00000000) T0 V2(b,x,w,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_b,s_x) == IEOR(s_w,s_a)) then
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) =  &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + 4.00000000d+00 & 
  * T0 & 
  * V2_(s_a, s_w, s_x)%array(i_a, i_w, i_x)

! Flop count 
flops = flops + psym(I_LENGTH, I_C, s_x) &
 * psym(I_LENGTH, I_C, s_w) &
 * psym(I_LENGTH, I_V, s_a) * 2.0d+00

end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_g_no0_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_g_no1_x0_type1_eri_v &
  (sb, ib, T0, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: T0
real(kind=8), intent(inout) :: V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_g_no1_x0_type1_eri_v &
  (sb, ib, T0, h2_i, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)

end subroutine g_if_sigma_ccvv_g_no1_x0_type1_eri_v



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_ccvv_g_no1_x0_type1_eri_v &
  (s_b, i_b, T0, V2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: T0
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_x, i_x, s_a, i_a
! S2(w,x,a,b) += (   -2.00000000) T0 V2(b,w,x,a) 
do s_w = 0, nir-1
do s_x = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_b,s_w) == IEOR(s_x,s_a)) then
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) =  &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  - 2.00000000d+00 & 
  * T0 & 
  * V2_(s_a, s_x, s_w)%array(i_a, i_x, i_w)

! Flop count 
flops = flops + psym(I_LENGTH, I_C, s_w) &
 * psym(I_LENGTH, I_C, s_x) &
 * psym(I_LENGTH, I_V, s_a) * 2.0d+00

end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_g_no1_x0_type1_eri_v

