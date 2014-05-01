#include <sci/icmr/fsrc/f_mr.fh>


!      ______                  __           
!     / ____/___   ____ ___   / /_ ____     
!    / /_   / _ \ / __ `__ \ / __// __ \ 
!   / __/  /  __// / / / / // /_ / /_/ /    
!  /_/     \___//_/ /_/ /_/ \__/ \____/  

!                                    Generated date : Sun Apr 20 10:26:23 2014



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ooov_no0_x0_type1_eri_c &
  (sx, ix, V2, W3, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sx, ix
real(kind=8), intent(inout) :: V2(*), W3(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sx, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sx

call set_symblock_Xcaa(sleft, W3, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_ccov_ooov_no0_x0_type1_eri_c &
  (sx, ix, h2_i, Xcaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_ccov_ooov_no0_x0_type1_eri_c



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

! -- Check1 is skipped 
!    -- is4RDM.first true
!    -- allRDM.first true
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ooov_no0_x0_type1_eri_c &
  (s_x, i_x, V2_, W3_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W3_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_w, i_w, s_a2, i_a2, s_a1, i_a1, s_a0, i_a0
! W3(x,w,a1,a0) += (    1.00000000) V2(x,a3,w,a2) D2(a3,a1,a2,a0) 
do s_a3 = 0, nir-1
do s_w = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_w) == IEOR(s_a1,s_a0) .and. & 
IEOR(s_x,s_a3) == IEOR(s_w,s_a2) .and. &
IEOR(s_a3,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- D2(a3,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_a0, i_a3, i_a2) =  &
  D2_(s_a0, s_a2, s_a1, s_a3)%array(i_a0, i_a2, i_a1, i_a3)
end do
end do
end do
end do
! Z2 <-- V2(x,a3,w,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a2, i_w) =  &
  V2_(s_a2, s_w, s_a3)%array(i_a2, i_w, i_a3)
end do
end do
end do

! Z3 <-- W3(x,w,a1,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0))

! W3(x,w,a1,a0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W3_(s_a0, s_a1, s_w)%array(i_a0, i_a1, i_w) = &
    W3_(s_a0, s_a1, s_w)%array(i_a0, i_a1, i_w) &
  + Z3_(i_a1, i_a0, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ooov_no0_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ooov_no0_x1_type1_eri_c &
  (sa, ia, sx, ix, T2, W3, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sx, ix
real(kind=8), intent(inout) :: T2(*), W3(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sx

call set_symblock_Xcaa(sleft, W3, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ooov_no0_x1_type1_eri_c &
  (sa, ia, sx, ix, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_ccov_ooov_no0_x1_type1_eri_c



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

! -- Check1 is skipped 
!    -- is4RDM.first false
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ooov_no0_x1_type1_eri_c &
  (s_a, i_a, s_x, i_x, T2_, W3_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W3_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_a1, i_a1, s_i, i_i, s_w, i_w
! S2(w,x,i,a) += (    1.00000000) T2(a0,a1,i,a) W3(x,w,a1,a0) 
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a0,s_a1) == IEOR(s_i,s_a) .and. &
IEOR(s_x,s_w) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- T2(a0,a1,i,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a1) =  &
  T2_(s_i, s_a1, s_a0)%array(i_i, i_a1, i_a0)
end do
end do
end do
! Z2 <-- W3(x,w,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a1, i_w) =  &
  W3_(s_a0, s_a1, s_w)%array(i_a0, i_a1, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_i, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ooov_no0_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ooov_no0_x0_type1_eri_o &
  (sa, ia, sa4, ia4, T2, W0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa4, ia4
real(kind=8), intent(inout) :: T2(*), W0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa4,sa)

call set_symblock_Xaa(sleft, W0, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_ccov_ooov_no0_x0_type1_eri_o &
  (sa, ia, sa4, ia4, av2_i, Xaa, d3, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ooov_no0_x0_type1_eri_o



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

! -- Check1 is skipped 
!    -- is4RDM.first false
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ooov_no0_x0_type1_eri_o &
  (s_a, i_a, s_a4, i_a4, T2_, W0_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a4, s_a4
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W0_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a2, i_a2, s_a0, i_a0, s_i, i_i, s_a3, i_a3
! W0(i,a3,a4,a) += (    1.00000000) T2(a1,a2,a0,a) D3(i,a3,a2,a4,a1,a0) 
do s_a1 = 0, nir-1
do s_a2 = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_a3 = 0, nir-1
if( &
IEOR(s_i,s_a3) == IEOR(s_a4,s_a) .and. & 
IEOR(s_a1,s_a2) == IEOR(s_a0,s_a) .and. &
IEOR(IEOR(s_i,s_a3),s_a2) == IEOR(IEOR(s_a4,s_a1),s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D3(i,a3,a2,a4,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a3, i_a2, i_a1, i_a0) =  &
  D3_(s_a0, s_a1, s_a4, s_a2, s_a3, s_i)%array(i_a0, i_a1, i_a4, i_a2, i_a3, i_i)
end do
end do
end do
end do
end do
! Z2 <-- T2(a1,a2,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1, i_a0) =  &
  T2_(s_a0, s_a2, s_a1)%array(i_a0, i_a2, i_a1)
end do
end do
end do

! Z3 <-- W0(i,a3,a4,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3),&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3))

! W0(i,a3,a4,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
W0_(s_a3, s_i)%array(i_a3, i_i) = &
    W0_(s_a3, s_i)%array(i_a3, i_i) &
  + Z3_(i_i, i_a3)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3) * &
                1 * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ooov_no0_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ooov_no0_x1_type1_eri_o &
  (sa, ia, sa4, ia4, V2, W0, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa4, ia4
real(kind=8), intent(inout) :: V2(*), W0(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa4, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa4,sa)

call set_symblock_Xaa(sleft, W0, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ooov_no0_x1_type1_eri_o &
  (sa, ia, sa4, ia4, h2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ooov_no0_x1_type1_eri_o



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

! -- Check1 is skipped 
!    -- is4RDM.first false
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ooov_no0_x1_type1_eri_o &
  (s_a, i_a, s_a4, i_a4, V2_, W0_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a4, s_a4
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W0_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a3, i_a3, s_i, i_i
! S2(w,x,i,a) += (    1.00000000) V2(a4,x,w,a3) W0(i,a3,a4,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a3 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a4,s_x) == IEOR(s_w,s_a3) .and. &
IEOR(s_i,s_a3) == IEOR(s_a4,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3) > 0) then

! Z1 <-- W0(i,a3,a4,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a3) =  &
  W0_(s_a3, s_i)%array(i_a3, i_i)
end do
end do
! Z2 <-- V2(a4,x,w,a3) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_x, i_w) =  &
  V2_(s_a3, s_w, s_x)%array(i_a3, i_w, i_x)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a3),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_i, i_x, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a3) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ooov_no0_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ooov_no1_x0_type1_eri_o &
  (sa, ia, sa3, ia3, T2, W1, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa3, ia3
real(kind=8), intent(inout) :: T2(*), W1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa3,sa)

call g_sigma_ccov_ooov_no1_x0_type1_eri_o &
  (sa, ia, sa3, ia3, av2_i, W1, d2, nir, nsym, psym, flops)

deallocate(av2_i)

end subroutine g_if_sigma_ccov_ooov_no1_x0_type1_eri_o



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

! -- Check1 is skipped 
!    -- is4RDM.first false
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second true
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ooov_no1_x0_type1_eri_o &
  (s_a, i_a, s_a3, i_a3, T2_, W1_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a3, s_a3
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W1_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_a1, i_a1, s_a2, i_a2
! W1(a3,a) += (    1.00000000) T2(a0,a1,a2,a) D2(a3,a1,a2,a0) 
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a2 = 0, nir-1
if( &
IEOR(s_a3,s_a) == 0 .and. & 
IEOR(s_a0,s_a1) == IEOR(s_a2,s_a) .and. &
IEOR(s_a3,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- T2(a0,a1,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1, i_a2) =  &
  T2_(s_a2, s_a1, s_a0)%array(i_a2, i_a1, i_a0)
end do
end do
end do
! Z2 <-- D2(a3,a1,a2,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a1, i_a2) =  &
  D2_(s_a0, s_a2, s_a1, s_a3)%array(i_a0, i_a2, i_a1, i_a3)
end do
end do
end do

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', 1,&
                     1,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     1,&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     1)

! W1(a3,a)  <-- Z3
W1_ = &
    W1_ &
  + Z3_

! Flop count
flops = flops + 1 * &
                1 * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ooov_no1_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ooov_no1_x1_type1_eri_o &
  (sa, ia, sa3, ia3, V2, W1, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa3, ia3
real(kind=8), intent(inout) :: V2(*), W1, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa3, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa3,sa)

call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ooov_no1_x1_type1_eri_o &
  (sa, ia, sa3, ia3, h2_i, W1, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)

end subroutine g_if_sigma_ccov_ooov_no1_x1_type1_eri_o



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

! -- Check1 is skipped 
!    -- is4RDM.first false
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ooov_no1_x1_type1_eri_o &
  (s_a, i_a, s_a3, i_a3, V2_, W1_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a3, s_a3
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W1_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_x, i_x, s_i, i_i
! S2(w,x,i,a) += (    1.00000000) V2(a3,w,x,i) W1(a3,a) 
do s_w = 0, nir-1
do s_x = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a3,s_w) == IEOR(s_x,s_i) .and. &
IEOR(s_a3,s_a) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) > 0) then

! Z1 <-- V2(a3,w,x,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_x, i_i) =  &
  V2_(s_i, s_x, s_w)%array(i_i, i_x, i_w)
end do
end do
end do
! Z2 <-- W1(a3,a) 
Z2_ =  &
  W1_

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     1,&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_w, i_x, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) * &
                1 * &
                1 * 2.0d+00

deallocate(Z1_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ooov_no1_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ooov_no2_x0_type1_eri_o &
  (sa, ia, sa3, ia3, T2, W2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa3, ia3
real(kind=8), intent(inout) :: T2(*), W2
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa3,sa)

call g_sigma_ccov_ooov_no2_x0_type1_eri_o &
  (sa, ia, sa3, ia3, av2_i, W2, d2, nir, nsym, psym, flops)

deallocate(av2_i)

end subroutine g_if_sigma_ccov_ooov_no2_x0_type1_eri_o



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

! -- Check1 is skipped 
!    -- is4RDM.first false
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second true
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ooov_no2_x0_type1_eri_o &
  (s_a, i_a, s_a3, i_a3, T2_, W2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a3, s_a3
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W2_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_a1, i_a1, s_a2, i_a2
! W2(a3,a) += (    1.00000000) T2(a0,a1,a2,a) D2(a3,a1,a2,a0) 
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a2 = 0, nir-1
if( &
IEOR(s_a3,s_a) == 0 .and. & 
IEOR(s_a0,s_a1) == IEOR(s_a2,s_a) .and. &
IEOR(s_a3,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- T2(a0,a1,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1, i_a2) =  &
  T2_(s_a2, s_a1, s_a0)%array(i_a2, i_a1, i_a0)
end do
end do
end do
! Z2 <-- D2(a3,a1,a2,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a1, i_a2) =  &
  D2_(s_a0, s_a2, s_a1, s_a3)%array(i_a0, i_a2, i_a1, i_a3)
end do
end do
end do

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', 1,&
                     1,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     1,&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     1)

! W2(a3,a)  <-- Z3
W2_ = &
    W2_ &
  + Z3_

! Flop count
flops = flops + 1 * &
                1 * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ooov_no2_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ooov_no2_x1_type1_eri_o &
  (sa, ia, sa3, ia3, V2, W2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa3, ia3
real(kind=8), intent(inout) :: V2(*), W2, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa3, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa3,sa)

call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ooov_no2_x1_type1_eri_o &
  (sa, ia, sa3, ia3, h2_i, W2, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)

end subroutine g_if_sigma_ccov_ooov_no2_x1_type1_eri_o



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

! -- Check1 is skipped 
!    -- is4RDM.first false
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ooov_no2_x1_type1_eri_o &
  (s_a, i_a, s_a3, i_a3, V2_, W2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a3, s_a3
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W2_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_i, i_i
! S2(w,x,i,a) += (   -2.00000000) V2(a3,x,w,i) W2(a3,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a3,s_x) == IEOR(s_w,s_i) .and. &
IEOR(s_a3,s_a) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0) then

! Z1 <-- V2(a3,x,w,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_w, i_i) =  &
  V2_(s_i, s_w, s_x)%array(i_i, i_w, i_x)
end do
end do
end do
! Z2 <-- W2(a3,a) 
Z2_ =  &
  W2_

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     1,&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_w, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                1 * &
                1 * 2.0d+00

deallocate(Z1_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ooov_no2_x1_type1_eri_o

