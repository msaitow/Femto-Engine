#include <sci/icmr/fsrc/f_mr.fh>


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

!                                    Generated date : Sun Apr 20 10:26:20 2014



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no0_x0_type0_noeri &
  (sb, ib, Fc0, T2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: Fc0
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no0_x0_type0_noeri &
  (sb, ib, Fc0, av2_i, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_covv_covv_no0_x0_type0_noeri



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
!    -- is4RDM.second true
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no0_x0_type0_noeri &
  (s_b, i_b, Fc0, T2_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Fc0
! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_a, i_a, s_i, i_i
! S2(w,i,a,b) += (    2.00000000) Fc0 T2(w,a0,a,b) D1(i,a0) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_a,s_b) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(w,a0,a,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  T2_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w)
end do
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00*Fc0, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no0_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no1_x0_type0_noeri &
  (sb, ib, Fc0, T2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: Fc0
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no1_x0_type0_noeri &
  (sb, ib, Fc0, av2_i, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_covv_covv_no1_x0_type0_noeri



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
!    -- is4RDM.second true
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no1_x0_type0_noeri &
  (s_b, i_b, Fc0, T2_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Fc0
! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_a, i_a, s_i, i_i
! S2(w,i,a,b) += (   -1.00000000) Fc0 T2(a0,w,a,b) D1(i,a0) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a,s_b) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,w,a,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  T2_(s_a, s_w, s_a0)%array(i_a, i_w, i_a0)
end do
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00*Fc0, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no1_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no2_x0_type0_noeri &
  (sb, ib, T2, W0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: T2(*), W0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xcav(sleft, W0, nir, nsym, psym) ! -> Xcav (allocate) 
call g_sigma_covv_covv_no2_x0_type0_noeri &
  (sb, ib, av2_i, Xcav, d1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no2_x0_type0_noeri



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
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no2_x0_type0_noeri &
  (s_b, i_b, T2_, W0_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a0, i_a0, s_a, i_a, s_i, i_i
! W0(c0,i,a,b) += (    1.00000000) T2(c0,a0,a,b) D1(i,a0) 
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_c0,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_c0,s_a0) == IEOR(s_a,s_b) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(c0,a0,a,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z1_(i_c0, i_a, i_a0) =  &
  T2_(s_a, s_a0, s_c0)%array(i_a, i_a0, i_c0)
end do
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- W0(c0,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_a))

! W0(c0,i,a,b)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W0_(s_a, s_i, s_c0)%array(i_a, i_i, i_c0) = &
    W0_(s_a, s_i, s_c0)%array(i_a, i_i, i_c0) &
  + Z3_(i_c0, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no2_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no2_x1_type0_noeri &
  (sb, ib, W0, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: W0(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sb

call set_symblock_Xcav(sleft, W0, nir, nsym, psym) ! -> Xcav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no2_x1_type0_noeri &
  (sb, ib, Xcav, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no2_x1_type0_noeri



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
subroutine g_sigma_covv_covv_no2_x1_type0_noeri &
  (s_b, i_b, W0_, S2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_i, i_i, s_a, i_a
! S2(w,i,a,b) += (   -2.00000000) Fc1(w,c0) W0(c0,i,a,b) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_c0) == 0 .and. &
IEOR(s_c0,s_i) == IEOR(s_a,s_b)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- W0(c0,i,a,b) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a, i_c0) =  &
  W0_(s_a, s_i, s_c0)%array(i_a, i_i, i_c0)
end do
end do
end do
! Z2 <-- Fc1(w,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  Fc1_(s_c0, s_w)%array(i_c0, i_w)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_i, i_a, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no2_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no3_x0_type0_noeri &
  (W1, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

real(kind=8), intent(inout) :: W1(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xaa(sleft, W1, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_covv_covv_no3_x0_type0_noeri &
  (Xaa, d2, fc1, nir, nsym, psym, flops)

deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no3_x0_type0_noeri



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
subroutine g_sigma_covv_covv_no3_x0_type0_noeri &
  (W1_, D2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W1_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_a2, i_a2, s_a1, i_a1
! W1(i,a0) += (    1.00000000) D2(i,a0,a2,a1) Fc1(a2,a1) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(s_i,s_a0) == 0 .and. & 
IEOR(s_i,s_a0) == IEOR(s_a2,s_a1) .and. &
IEOR(s_a2,s_a1) == 0) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(i,a0,a2,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a2, i_a1) =  &
  D2_(s_a1, s_a2, s_a0, s_i)%array(i_a1, i_a2, i_a0, i_i)
end do
end do
end do
end do
! Z2 <-- Fc1(a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1) =  &
  Fc1_(s_a1, s_a2)%array(i_a1, i_a2)
end do
end do

! Z3 <-- W1(i,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W1(i,a0)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W1_(s_a0, s_i)%array(i_a0, i_i) = &
    W1_(s_a0, s_i)%array(i_a0, i_i) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no3_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no3_x1_type0_noeri &
  (sb, ib, T2, W1, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: T2(*), W1(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xaa(sleft, W1, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no3_x1_type0_noeri &
  (sb, ib, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no3_x1_type0_noeri



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
subroutine g_sigma_covv_covv_no3_x1_type0_noeri &
  (s_b, i_b, T2_, W1_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W1_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_a, i_a, s_i, i_i
! S2(w,i,a,b) += (    2.00000000) T2(w,a0,a,b) W1(i,a0) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_a,s_b) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(w,a0,a,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  T2_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w)
end do
end do
end do
! Z2 <-- W1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  W1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no3_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no4_x0_type0_noeri &
  (sb, ib, T2, W2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: T2(*), W2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xcav(sleft, W2, nir, nsym, psym) ! -> Xcav (allocate) 
call g_sigma_covv_covv_no4_x0_type0_noeri &
  (sb, ib, av2_i, Xcav, d1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no4_x0_type0_noeri



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
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no4_x0_type0_noeri &
  (s_b, i_b, T2_, W2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_c0, i_c0, s_a, i_a, s_i, i_i
! W2(c0,i,b,a) += (    1.00000000) T2(a0,c0,a,b) D1(i,a0) 
do s_a0 = 0, nir-1
do s_c0 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_c0,s_i) == IEOR(s_b,s_a) .and. & 
IEOR(s_a0,s_c0) == IEOR(s_a,s_b) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,c0,a,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z1_(i_c0, i_a, i_a0) =  &
  T2_(s_a, s_c0, s_a0)%array(i_a, i_c0, i_a0)
end do
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- W2(c0,i,b,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_a))

! W2(c0,i,b,a)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W2_(s_a, s_i, s_c0)%array(i_a, i_i, i_c0) = &
    W2_(s_a, s_i, s_c0)%array(i_a, i_i, i_c0) &
  + Z3_(i_c0, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no4_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no4_x1_type0_noeri &
  (sb, ib, W2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: W2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sb

call set_symblock_Xcav(sleft, W2, nir, nsym, psym) ! -> Xcav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no4_x1_type0_noeri &
  (sb, ib, Xcav, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no4_x1_type0_noeri



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
subroutine g_sigma_covv_covv_no4_x1_type0_noeri &
  (s_b, i_b, W2_, S2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_i, i_i, s_a, i_a
! S2(w,i,a,b) += (    1.00000000) Fc1(w,c0) W2(c0,i,b,a) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_c0) == 0 .and. &
IEOR(s_c0,s_i) == IEOR(s_b,s_a)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- W2(c0,i,b,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a, i_c0) =  &
  W2_(s_a, s_i, s_c0)%array(i_a, i_i, i_c0)
end do
end do
end do
! Z2 <-- Fc1(w,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  Fc1_(s_c0, s_w)%array(i_c0, i_w)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_i, i_a, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no4_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no5_x0_type0_noeri &
  (W3, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

real(kind=8), intent(inout) :: W3(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xaa(sleft, W3, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_covv_covv_no5_x0_type0_noeri &
  (Xaa, d2, fc1, nir, nsym, psym, flops)

deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no5_x0_type0_noeri



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
subroutine g_sigma_covv_covv_no5_x0_type0_noeri &
  (W3_, D2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W3_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_a2, i_a2, s_a1, i_a1
! W3(i,a0) += (    1.00000000) D2(i,a0,a2,a1) Fc1(a2,a1) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(s_i,s_a0) == 0 .and. & 
IEOR(s_i,s_a0) == IEOR(s_a2,s_a1) .and. &
IEOR(s_a2,s_a1) == 0) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(i,a0,a2,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a2, i_a1) =  &
  D2_(s_a1, s_a2, s_a0, s_i)%array(i_a1, i_a2, i_a0, i_i)
end do
end do
end do
end do
! Z2 <-- Fc1(a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1) =  &
  Fc1_(s_a1, s_a2)%array(i_a1, i_a2)
end do
end do

! Z3 <-- W3(i,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W3(i,a0)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W3_(s_a0, s_i)%array(i_a0, i_i) = &
    W3_(s_a0, s_i)%array(i_a0, i_i) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no5_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no5_x1_type0_noeri &
  (sb, ib, T2, W3, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: T2(*), W3(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xaa(sleft, W3, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no5_x1_type0_noeri &
  (sb, ib, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no5_x1_type0_noeri



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
subroutine g_sigma_covv_covv_no5_x1_type0_noeri &
  (s_b, i_b, T2_, W3_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W3_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_a, i_a, s_i, i_i
! S2(w,i,a,b) += (   -1.00000000) T2(a0,w,a,b) W3(i,a0) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a,s_b) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,w,a,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  T2_(s_a, s_w, s_a0)%array(i_a, i_w, i_a0)
end do
end do
end do
! Z2 <-- W3(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  W3_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no5_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no6_x0_type0_noeri &
  (sa, ia, T2, W4, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), W4(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcav(sleft, W4, nir, nsym, psym) ! -> Xcav (allocate) 
call g_sigma_covv_covv_no6_x0_type0_noeri &
  (sa, ia, av2_i, Xcav, d1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no6_x0_type0_noeri



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
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no6_x0_type0_noeri &
  (s_a, i_a, T2_, W4_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W4_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_v0, i_v0, s_i, i_i
! W4(w,i,v0,a) += (    1.00000000) T2(w,a0,v0,a) D1(i,a0) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_v0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_v0,s_a) .and. & 
IEOR(s_w,s_a0) == IEOR(s_v0,s_a) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(w,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_v0, i_a0) =  &
  T2_(s_v0, s_a0, s_w)%array(i_v0, i_a0, i_w)
end do
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- W4(w,i,v0,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0))

! W4(w,i,v0,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
W4_(s_v0, s_i, s_w)%array(i_v0, i_i, i_w) = &
    W4_(s_v0, s_i, s_w)%array(i_v0, i_i, i_w) &
  + Z3_(i_w, i_v0, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no6_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no6_x1_type0_noeri &
  (sa, ia, sb, ib, W4, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sb, ib
real(kind=8), intent(inout) :: W4(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa

call set_symblock_Xcav(sleft, W4, nir, nsym, psym) ! -> Xcav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no6_x1_type0_noeri &
  (sa, ia, sb, ib, Xcav, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no6_x1_type0_noeri



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
subroutine g_sigma_covv_covv_no6_x1_type0_noeri &
  (s_a, i_a, s_b, i_b, W4_, S2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W4_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_w, i_w, s_i, i_i
! S2(w,i,a,b) += (   -1.00000000) Fc1(v0,b) W4(w,i,v0,a) 
do s_v0 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_v0,s_b) == 0 .and. &
IEOR(s_w,s_i) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W4(w,i,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i, i_v0) =  &
  W4_(s_v0, s_i, s_w)%array(i_v0, i_i, i_w)
end do
end do
end do
! Z2 <-- Fc1(v0,b) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0) =  &
  Fc1_(s_b, s_v0)%array(i_b, i_v0)
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_V, s_v0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                1 * &
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no6_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no7_x0_type0_noeri &
  (sb, ib, T2, W5, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: T2(*), W5(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xcav(sleft, W5, nir, nsym, psym) ! -> Xcav (allocate) 
call g_sigma_covv_covv_no7_x0_type0_noeri &
  (sb, ib, av2_i, Xcav, d1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no7_x0_type0_noeri



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
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no7_x0_type0_noeri &
  (s_b, i_b, T2_, W5_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W5_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_v0, i_v0, s_i, i_i
! W5(w,i,v0,b) += (    1.00000000) T2(w,a0,v0,b) D1(i,a0) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_v0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_v0,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_v0,s_b) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(w,a0,v0,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_v0, i_a0) =  &
  T2_(s_v0, s_a0, s_w)%array(i_v0, i_a0, i_w)
end do
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- W5(w,i,v0,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0))

! W5(w,i,v0,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
W5_(s_v0, s_i, s_w)%array(i_v0, i_i, i_w) = &
    W5_(s_v0, s_i, s_w)%array(i_v0, i_i, i_w) &
  + Z3_(i_w, i_v0, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no7_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no7_x1_type0_noeri &
  (sb, ib, W5, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: W5(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sb

call set_symblock_Xcav(sleft, W5, nir, nsym, psym) ! -> Xcav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no7_x1_type0_noeri &
  (sb, ib, Xcav, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no7_x1_type0_noeri



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
subroutine g_sigma_covv_covv_no7_x1_type0_noeri &
  (s_b, i_b, W5_, S2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W5_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_a, i_a, s_w, i_w, s_i, i_i
! S2(w,i,a,b) += (    2.00000000) Fc1(v0,a) W5(w,i,v0,b) 
do s_v0 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_v0,s_a) == 0 .and. &
IEOR(s_w,s_i) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- Fc1(v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_v0) =  &
  Fc1_(s_a, s_v0)%array(i_a, i_v0)
end do
end do
! Z2 <-- W5(w,i,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_w, i_i) =  &
  W5_(s_v0, s_i, s_w)%array(i_v0, i_i, i_w)
end do
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_V, s_v0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_a, i_w, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no7_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no8_x0_type0_noeri &
  (sv0, iv0, T2, W6, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), W6(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xcav(sleft, W6, nir, nsym, psym) ! -> Xcav (allocate) 
call g_sigma_covv_covv_no8_x0_type0_noeri &
  (sv0, iv0, av2_i, Xcav, d1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no8_x0_type0_noeri



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
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no8_x0_type0_noeri &
  (s_v0, i_v0, T2_, W6_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W6_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_a, i_a, s_i, i_i
! W6(w,i,a,v0) += (    1.00000000) T2(w,a0,a,v0) D1(i,a0) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_v0) .and. & 
IEOR(s_w,s_a0) == IEOR(s_a,s_v0) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(w,a0,a,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  T2_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w)
end do
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- W6(w,i,a,v0) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! W6(w,i,a,v0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W6_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    W6_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no8_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no8_x1_type0_noeri &
  (sb, ib, sv0, iv0, W6, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: W6(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sv0

call set_symblock_Xcav(sleft, W6, nir, nsym, psym) ! -> Xcav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no8_x1_type0_noeri &
  (sb, ib, sv0, iv0, Xcav, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no8_x1_type0_noeri



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
subroutine g_sigma_covv_covv_no8_x1_type0_noeri &
  (s_b, i_b, s_v0, i_v0, W6_, S2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W6_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_i, i_i, s_a, i_a
! S2(w,i,a,b) += (    2.00000000) Fc1(v0,b) W6(w,i,a,v0) 
do s_w = 0, nir-1
do s_i = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_v0,s_b) == 0 .and. &
IEOR(s_w,s_i) == IEOR(s_a,s_v0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) > 0) then

! Z1 <-- W6(w,i,a,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a)))
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i, i_a) =  &
  W6_(s_a, s_i, s_w)%array(i_a, i_i, i_w)
end do
end do
end do
! Z2 <-- Fc1(v0,b) 
Z2_ =  &
  Fc1_(s_b, s_v0)%array(i_b, i_v0)

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     1,&
                     1,&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_i, i_a)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) * &
                1 * &
                1 * 2.0d+00

deallocate(Z1_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no8_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no9_x0_type0_noeri &
  (sb, ib, T2, W7, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: T2(*), W7(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xcav(sleft, W7, nir, nsym, psym) ! -> Xcav (allocate) 
call g_sigma_covv_covv_no9_x0_type0_noeri &
  (sb, ib, av2_i, Xcav, d1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no9_x0_type0_noeri



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
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no9_x0_type0_noeri &
  (s_b, i_b, T2_, W7_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W7_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_v0, i_v0, s_i, i_i
! W7(w,i,b,v0) += (    1.00000000) T2(a0,w,v0,b) D1(i,a0) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_v0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_b,s_v0) .and. & 
IEOR(s_a0,s_w) == IEOR(s_v0,s_b) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,w,v0,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_v0, i_a0) =  &
  T2_(s_v0, s_w, s_a0)%array(i_v0, i_w, i_a0)
end do
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- W7(w,i,b,v0) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0))

! W7(w,i,b,v0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
W7_(s_v0, s_i, s_w)%array(i_v0, i_i, i_w) = &
    W7_(s_v0, s_i, s_w)%array(i_v0, i_i, i_w) &
  + Z3_(i_w, i_v0, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no9_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no9_x1_type0_noeri &
  (sb, ib, W7, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: W7(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sb

call set_symblock_Xcav(sleft, W7, nir, nsym, psym) ! -> Xcav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no9_x1_type0_noeri &
  (sb, ib, Xcav, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no9_x1_type0_noeri



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
subroutine g_sigma_covv_covv_no9_x1_type0_noeri &
  (s_b, i_b, W7_, S2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W7_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_a, i_a, s_w, i_w, s_i, i_i
! S2(w,i,a,b) += (   -1.00000000) Fc1(v0,a) W7(w,i,b,v0) 
do s_v0 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_v0,s_a) == 0 .and. &
IEOR(s_w,s_i) == IEOR(s_b,s_v0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- Fc1(v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_v0) =  &
  Fc1_(s_a, s_v0)%array(i_a, i_v0)
end do
end do
! Z2 <-- W7(w,i,b,v0) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_w, i_i) =  &
  W7_(s_v0, s_i, s_w)%array(i_v0, i_i, i_w)
end do
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_V, s_v0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_a, i_w, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no9_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no0_x0_type1_eri_c &
  (sw, iw, V2, W8, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), W8(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_Xcaa(sleft, W8, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_covv_covv_no0_x0_type1_eri_c &
  (sw, iw, h2_i, Xcaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_covv_covv_no0_x0_type1_eri_c



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
subroutine g_sigma_covv_covv_no0_x0_type1_eri_c &
  (s_w, i_w, V2_, W8_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W8_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_c0, i_c0, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W8(w,c0,i,a0) += (    1.00000000) V2(w,a2,c0,a1) D2(i,a0,a2,a1) 
do s_a2 = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_c0) == IEOR(s_i,s_a0) .and. & 
IEOR(s_w,s_a2) == IEOR(s_c0,s_a1) .and. &
IEOR(s_i,s_a0) == IEOR(s_a2,s_a1)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(i,a0,a2,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a2, i_a1) =  &
  D2_(s_a1, s_a2, s_a0, s_i)%array(i_a1, i_a2, i_a0, i_i)
end do
end do
end do
end do
! Z2 <-- V2(w,a2,c0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1, i_c0) =  &
  V2_(s_a1, s_c0, s_a2)%array(i_a1, i_c0, i_a2)
end do
end do
end do

! Z3 <-- W8(w,c0,i,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W8(w,c0,i,a0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W8_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0) = &
    W8_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0) &
  + Z3_(i_i, i_a0, i_c0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no0_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no0_x1_type1_eri_c &
  (sb, ib, sw, iw, T2, W8, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: T2(*), W8(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xcaa(sleft, W8, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no0_x1_type1_eri_c &
  (sb, ib, sw, iw, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_covv_covv_no0_x1_type1_eri_c



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
subroutine g_sigma_covv_covv_no0_x1_type1_eri_c &
  (s_b, i_b, s_w, i_w, T2_, W8_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W8_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a0, i_a0, s_a, i_a, s_i, i_i
! S2(w,i,a,b) += (    1.00000000) T2(c0,a0,a,b) W8(w,c0,i,a0) 
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_c0,s_a0) == IEOR(s_a,s_b) .and. &
IEOR(s_w,s_c0) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(c0,a0,a,b) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_c0, i_a0) =  &
  T2_(s_a, s_a0, s_c0)%array(i_a, i_a0, i_c0)
end do
end do
end do
! Z2 <-- W8(w,c0,i,a0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0, i_i) =  &
  W8_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0)
end do
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_a, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no0_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no1_x0_type1_eri_c &
  (sw, iw, V2, W9, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), W9(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_Xcaa(sleft, W9, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_covv_covv_no1_x0_type1_eri_c &
  (sw, iw, h2_i, Xcaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_covv_covv_no1_x0_type1_eri_c



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
subroutine g_sigma_covv_covv_no1_x0_type1_eri_c &
  (s_w, i_w, V2_, W9_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W9_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a2, i_a2, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W9(w,c0,i,a0) += (    1.00000000) V2(w,c0,a2,a1) D2(i,a0,a2,a1) 
do s_c0 = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_c0) == IEOR(s_i,s_a0) .and. & 
IEOR(s_w,s_c0) == IEOR(s_a2,s_a1) .and. &
IEOR(s_i,s_a0) == IEOR(s_a2,s_a1)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(i,a0,a2,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a2, i_a1) =  &
  D2_(s_a1, s_a2, s_a0, s_i)%array(i_a1, i_a2, i_a0, i_i)
end do
end do
end do
end do
! Z2 <-- V2(w,c0,a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1, i_c0) =  &
  V2_(s_a1, s_a2, s_c0)%array(i_a1, i_a2, i_c0)
end do
end do
end do

! Z3 <-- W9(w,c0,i,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W9(w,c0,i,a0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W9_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0) = &
    W9_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0) &
  + Z3_(i_i, i_a0, i_c0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no1_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no1_x1_type1_eri_c &
  (sb, ib, sw, iw, T2, W9, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: T2(*), W9(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xcaa(sleft, W9, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no1_x1_type1_eri_c &
  (sb, ib, sw, iw, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_covv_covv_no1_x1_type1_eri_c



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
subroutine g_sigma_covv_covv_no1_x1_type1_eri_c &
  (s_b, i_b, s_w, i_w, T2_, W9_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W9_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a0, i_a0, s_a, i_a, s_i, i_i
! S2(w,i,a,b) += (   -2.00000000) T2(c0,a0,a,b) W9(w,c0,i,a0) 
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_c0,s_a0) == IEOR(s_a,s_b) .and. &
IEOR(s_w,s_c0) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(c0,a0,a,b) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_c0, i_a0) =  &
  T2_(s_a, s_a0, s_c0)%array(i_a, i_a0, i_c0)
end do
end do
end do
! Z2 <-- W9(w,c0,i,a0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0, i_i) =  &
  W9_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0)
end do
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_a, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no1_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no2_x0_type1_eri_c &
  (sw, iw, V2, W11, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), W11(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_Xcaa(sleft, W11, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_covv_covv_no2_x0_type1_eri_c &
  (sw, iw, h2_i, Xcaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_covv_covv_no2_x0_type1_eri_c



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
subroutine g_sigma_covv_covv_no2_x0_type1_eri_c &
  (s_w, i_w, V2_, W11_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W11_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_c0, i_c0, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W11(w,c0,i,a0) += (    1.00000000) V2(w,a2,c0,a1) D2(i,a1,a2,a0) 
do s_a2 = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_c0) == IEOR(s_i,s_a0) .and. & 
IEOR(s_w,s_a2) == IEOR(s_c0,s_a1) .and. &
IEOR(s_i,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- D2(i,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a1, i_a2) =  &
  D2_(s_a0, s_a2, s_a1, s_i)%array(i_a0, i_a2, i_a1, i_i)
end do
end do
end do
end do
! Z2 <-- V2(w,a2,c0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a2, i_c0) =  &
  V2_(s_a1, s_c0, s_a2)%array(i_a1, i_c0, i_a2)
end do
end do
end do

! Z3 <-- W11(w,c0,i,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W11(w,c0,i,a0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W11_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0) = &
    W11_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0) &
  + Z3_(i_i, i_a0, i_c0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no2_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no2_x1_type1_eri_c &
  (sb, ib, sw, iw, T2, W11, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: T2(*), W11(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xcaa(sleft, W11, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no2_x1_type1_eri_c &
  (sb, ib, sw, iw, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_covv_covv_no2_x1_type1_eri_c



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
subroutine g_sigma_covv_covv_no2_x1_type1_eri_c &
  (s_b, i_b, s_w, i_w, T2_, W11_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W11_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_c0, i_c0, s_a, i_a, s_i, i_i
! S2(w,i,a,b) += (    1.00000000) T2(a0,c0,a,b) W11(w,c0,i,a0) 
do s_a0 = 0, nir-1
do s_c0 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_c0) == IEOR(s_a,s_b) .and. &
IEOR(s_w,s_c0) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(a0,c0,a,b) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0, i_c0) =  &
  T2_(s_a, s_c0, s_a0)%array(i_a, i_c0, i_a0)
end do
end do
end do
! Z2 <-- W11(w,c0,i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_c0, i_i) =  &
  W11_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0)
end do
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_a, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no2_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no3_x0_type1_eri_c &
  (sw, iw, V2, W12, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), W12(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_Xcaa(sleft, W12, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_covv_covv_no3_x0_type1_eri_c &
  (sw, iw, h2_i, Xcaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_covv_covv_no3_x0_type1_eri_c



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
subroutine g_sigma_covv_covv_no3_x0_type1_eri_c &
  (s_w, i_w, V2_, W12_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W12_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a2, i_a2, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W12(w,c0,i,a0) += (    1.00000000) V2(w,c0,a2,a1) D2(i,a0,a2,a1) 
do s_c0 = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_c0) == IEOR(s_i,s_a0) .and. & 
IEOR(s_w,s_c0) == IEOR(s_a2,s_a1) .and. &
IEOR(s_i,s_a0) == IEOR(s_a2,s_a1)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(i,a0,a2,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a2, i_a1) =  &
  D2_(s_a1, s_a2, s_a0, s_i)%array(i_a1, i_a2, i_a0, i_i)
end do
end do
end do
end do
! Z2 <-- V2(w,c0,a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1, i_c0) =  &
  V2_(s_a1, s_a2, s_c0)%array(i_a1, i_a2, i_c0)
end do
end do
end do

! Z3 <-- W12(w,c0,i,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W12(w,c0,i,a0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W12_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0) = &
    W12_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0) &
  + Z3_(i_i, i_a0, i_c0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no3_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no3_x1_type1_eri_c &
  (sb, ib, sw, iw, T2, W12, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: T2(*), W12(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xcaa(sleft, W12, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no3_x1_type1_eri_c &
  (sb, ib, sw, iw, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_covv_covv_no3_x1_type1_eri_c



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
subroutine g_sigma_covv_covv_no3_x1_type1_eri_c &
  (s_b, i_b, s_w, i_w, T2_, W12_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W12_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_c0, i_c0, s_a, i_a, s_i, i_i
! S2(w,i,a,b) += (    1.00000000) T2(a0,c0,a,b) W12(w,c0,i,a0) 
do s_a0 = 0, nir-1
do s_c0 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_c0) == IEOR(s_a,s_b) .and. &
IEOR(s_w,s_c0) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(a0,c0,a,b) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0, i_c0) =  &
  T2_(s_a, s_c0, s_a0)%array(i_a, i_c0, i_a0)
end do
end do
end do
! Z2 <-- W12(w,c0,i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_c0, i_i) =  &
  W12_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0)
end do
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_a, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no3_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no4_x0_type1_eri_c &
  (sb, ib, sw, iw, T2, V2, W14, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: T2(*), V2(*), W14(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sw,sb)

call set_symblock_Xav(sleft, W14, nir, nsym, psym) ! -> Xav (allocate) 
call g_sigma_covv_covv_no4_x0_type1_eri_c &
  (sb, ib, sw, iw, av2_i, h2_i, Xav, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xav)

end subroutine g_if_sigma_covv_covv_no4_x0_type1_eri_c



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
subroutine g_sigma_covv_covv_no4_x0_type1_eri_c &
  (s_b, i_b, s_w, i_w, T2_, V2_, W14_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W14_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_v0, i_v0, s_c0, i_c0, s_a0, i_a0
! W14(w,a0,b,a) += (    1.00000000) V2(w,a,v0,c0) T2(c0,a0,v0,b) 
do s_a = 0, nir-1
do s_v0 = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_b,s_a) .and. & 
IEOR(s_w,s_a) == IEOR(s_v0,s_c0) .and. &
IEOR(s_c0,s_a0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(w,a,v0,c0) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_v0, i_c0) =  &
  V2_(s_c0, s_v0, s_a)%array(i_c0, i_v0, i_a)
end do
end do
end do
! Z2 <-- T2(c0,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_c0, i_a0) =  &
  T2_(s_v0, s_a0, s_c0)%array(i_v0, i_a0, i_c0)
end do
end do
end do

! Z3 <-- W14(w,a0,b,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W14(w,a0,b,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W14_(s_a, s_a0)%array(i_a, i_a0) = &
    W14_(s_a, s_a0)%array(i_a, i_a0) &
  + Z3_(i_a, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no4_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no4_x1_type1_eri_c &
  (sb, ib, sw, iw, W14, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: W14(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sw,sb)

call set_symblock_Xav(sleft, W14, nir, nsym, psym) ! -> Xav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no4_x1_type1_eri_c &
  (sb, ib, sw, iw, Xav, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xav)

end subroutine g_if_sigma_covv_covv_no4_x1_type1_eri_c



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
!    -- is4RDM.second true
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no4_x1_type1_eri_c &
  (s_b, i_b, s_w, i_w, W14_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W14_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_a, i_a
! S2(w,i,a,b) += (    4.00000000) D1(i,a0) W14(w,a0,b,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_b,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W14(w,a0,b,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0) =  &
  W14_(s_a, s_a0)%array(i_a, i_a0)
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_a, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no4_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no5_x0_type1_eri_c &
  (sb, ib, sw, iw, T2, V2, W20, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: T2(*), V2(*), W20(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sw,sb)

call set_symblock_Xav(sleft, W20, nir, nsym, psym) ! -> Xav (allocate) 
call g_sigma_covv_covv_no5_x0_type1_eri_c &
  (sb, ib, sw, iw, av2_i, h2_i, Xav, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xav)

end subroutine g_if_sigma_covv_covv_no5_x0_type1_eri_c



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
subroutine g_sigma_covv_covv_no5_x0_type1_eri_c &
  (s_b, i_b, s_w, i_w, T2_, V2_, W20_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W20_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_v0, i_v0, s_a, i_a, s_a0, i_a0
! W20(w,a0,b,a) += (    1.00000000) V2(w,c0,v0,a) T2(c0,a0,v0,b) 
do s_c0 = 0, nir-1
do s_v0 = 0, nir-1
do s_a = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_b,s_a) .and. & 
IEOR(s_w,s_c0) == IEOR(s_v0,s_a) .and. &
IEOR(s_c0,s_a0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- V2(w,c0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_c0, i_v0) =  &
  V2_(s_a, s_v0, s_c0)%array(i_a, i_v0, i_c0)
end do
end do
end do
! Z2 <-- T2(c0,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_v0, i_a0) =  &
  T2_(s_v0, s_a0, s_c0)%array(i_v0, i_a0, i_c0)
end do
end do
end do

! Z3 <-- W20(w,a0,b,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W20(w,a0,b,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W20_(s_a, s_a0)%array(i_a, i_a0) = &
    W20_(s_a, s_a0)%array(i_a, i_a0) &
  + Z3_(i_a, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no5_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no5_x1_type1_eri_c &
  (sb, ib, sw, iw, W20, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: W20(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sw,sb)

call set_symblock_Xav(sleft, W20, nir, nsym, psym) ! -> Xav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no5_x1_type1_eri_c &
  (sb, ib, sw, iw, Xav, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xav)

end subroutine g_if_sigma_covv_covv_no5_x1_type1_eri_c



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
!    -- is4RDM.second true
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no5_x1_type1_eri_c &
  (s_b, i_b, s_w, i_w, W20_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W20_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_a, i_a
! S2(w,i,a,b) += (   -2.00000000) D1(i,a0) W20(w,a0,b,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_b,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W20(w,a0,b,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0) =  &
  W20_(s_a, s_a0)%array(i_a, i_a0)
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_a, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no5_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no6_x0_type1_eri_c &
  (sb, ib, sw, iw, T2, V2, W22, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: T2(*), V2(*), W22(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sw,sb)

call set_symblock_Xav(sleft, W22, nir, nsym, psym) ! -> Xav (allocate) 
call g_sigma_covv_covv_no6_x0_type1_eri_c &
  (sb, ib, sw, iw, av2_i, h2_i, Xav, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xav)

end subroutine g_if_sigma_covv_covv_no6_x0_type1_eri_c



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
subroutine g_sigma_covv_covv_no6_x0_type1_eri_c &
  (s_b, i_b, s_w, i_w, T2_, V2_, W22_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W22_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_v0, i_v0, s_c0, i_c0, s_a0, i_a0
! W22(w,a0,b,a) += (    1.00000000) V2(w,a,v0,c0) T2(a0,c0,v0,b) 
do s_a = 0, nir-1
do s_v0 = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_b,s_a) .and. & 
IEOR(s_w,s_a) == IEOR(s_v0,s_c0) .and. &
IEOR(s_a0,s_c0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(w,a,v0,c0) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_v0, i_c0) =  &
  V2_(s_c0, s_v0, s_a)%array(i_c0, i_v0, i_a)
end do
end do
end do
! Z2 <-- T2(a0,c0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_c0, i_a0) =  &
  T2_(s_v0, s_c0, s_a0)%array(i_v0, i_c0, i_a0)
end do
end do
end do

! Z3 <-- W22(w,a0,b,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W22(w,a0,b,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W22_(s_a, s_a0)%array(i_a, i_a0) = &
    W22_(s_a, s_a0)%array(i_a, i_a0) &
  + Z3_(i_a, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no6_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no6_x1_type1_eri_c &
  (sb, ib, sw, iw, W22, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: W22(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sw,sb)

call set_symblock_Xav(sleft, W22, nir, nsym, psym) ! -> Xav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no6_x1_type1_eri_c &
  (sb, ib, sw, iw, Xav, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xav)

end subroutine g_if_sigma_covv_covv_no6_x1_type1_eri_c



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
!    -- is4RDM.second true
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no6_x1_type1_eri_c &
  (s_b, i_b, s_w, i_w, W22_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W22_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_a, i_a
! S2(w,i,a,b) += (   -2.00000000) D1(i,a0) W22(w,a0,b,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_b,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W22(w,a0,b,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0) =  &
  W22_(s_a, s_a0)%array(i_a, i_a0)
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_a, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no6_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no7_x0_type1_eri_c &
  (sb, ib, sw, iw, T2, V2, W28, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: T2(*), V2(*), W28(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sw,sb)

call set_symblock_Xav(sleft, W28, nir, nsym, psym) ! -> Xav (allocate) 
call g_sigma_covv_covv_no7_x0_type1_eri_c &
  (sb, ib, sw, iw, av2_i, h2_i, Xav, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xav)

end subroutine g_if_sigma_covv_covv_no7_x0_type1_eri_c



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
subroutine g_sigma_covv_covv_no7_x0_type1_eri_c &
  (s_b, i_b, s_w, i_w, T2_, V2_, W28_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W28_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_v0, i_v0, s_a, i_a, s_a0, i_a0
! W28(w,a0,b,a) += (    1.00000000) V2(w,c0,v0,a) T2(a0,c0,v0,b) 
do s_c0 = 0, nir-1
do s_v0 = 0, nir-1
do s_a = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_b,s_a) .and. & 
IEOR(s_w,s_c0) == IEOR(s_v0,s_a) .and. &
IEOR(s_a0,s_c0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- V2(w,c0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_c0, i_v0) =  &
  V2_(s_a, s_v0, s_c0)%array(i_a, i_v0, i_c0)
end do
end do
end do
! Z2 <-- T2(a0,c0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_v0, i_a0) =  &
  T2_(s_v0, s_c0, s_a0)%array(i_v0, i_c0, i_a0)
end do
end do
end do

! Z3 <-- W28(w,a0,b,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W28(w,a0,b,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W28_(s_a, s_a0)%array(i_a, i_a0) = &
    W28_(s_a, s_a0)%array(i_a, i_a0) &
  + Z3_(i_a, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no7_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no7_x1_type1_eri_c &
  (sb, ib, sw, iw, W28, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: W28(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sw,sb)

call set_symblock_Xav(sleft, W28, nir, nsym, psym) ! -> Xav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no7_x1_type1_eri_c &
  (sb, ib, sw, iw, Xav, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xav)

end subroutine g_if_sigma_covv_covv_no7_x1_type1_eri_c



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
!    -- is4RDM.second true
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no7_x1_type1_eri_c &
  (s_b, i_b, s_w, i_w, W28_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W28_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_a, i_a
! S2(w,i,a,b) += (    1.00000000) D1(i,a0) W28(w,a0,b,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_b,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W28(w,a0,b,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0) =  &
  W28_(s_a, s_a0)%array(i_a, i_a0)
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_a, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no7_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no0_x0_type1_eri_o &
  (sa4, ia4, V2, W10, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa4, ia4
real(kind=8), intent(inout) :: V2(*), W10(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa4, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_Xaa(sleft, W10, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_covv_covv_no0_x0_type1_eri_o &
  (sa4, ia4, h2_i, Xaa, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no0_x0_type1_eri_o



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
subroutine g_sigma_covv_covv_no0_x0_type1_eri_o &
  (s_a4, i_a4, V2_, W10_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a4, s_a4
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W10_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a3, i_a3, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W10(i,a0) += (    1.00000000) V2(a4,a2,a3,a1) D3(i,a0,a4,a2,a3,a1) 
do s_a2 = 0, nir-1
do s_a3 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == 0 .and. & 
IEOR(s_a4,s_a2) == IEOR(s_a3,s_a1) .and. &
IEOR(IEOR(s_i,s_a0),s_a4) == IEOR(IEOR(s_a2,s_a3),s_a1)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D3(i,a0,a4,a2,a3,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a2, i_a3, i_a1) =  &
  D3_(s_a1, s_a3, s_a2, s_a4, s_a0, s_i)%array(i_a1, i_a3, i_a2, i_a4, i_a0, i_i)
end do
end do
end do
end do
end do
! Z2 <-- V2(a4,a2,a3,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a3, i_a1) =  &
  V2_(s_a1, s_a3, s_a2)%array(i_a1, i_a3, i_a2)
end do
end do
end do

! Z3 <-- W10(i,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W10(i,a0)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W10_(s_a0, s_i)%array(i_a0, i_i) = &
    W10_(s_a0, s_i)%array(i_a0, i_i) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no0_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no1_x0_type1_eri_o &
  (sa4, ia4, V2, W13, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa4, ia4
real(kind=8), intent(inout) :: V2(*), W13(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa4, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_Xaa(sleft, W13, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_covv_covv_no1_x0_type1_eri_o &
  (sa4, ia4, h2_i, Xaa, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no1_x0_type1_eri_o



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
subroutine g_sigma_covv_covv_no1_x0_type1_eri_o &
  (s_a4, i_a4, V2_, W13_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a4, s_a4
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W13_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a3, i_a3, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W13(i,a0) += (    1.00000000) V2(a4,a2,a3,a1) D3(i,a0,a4,a2,a3,a1) 
do s_a2 = 0, nir-1
do s_a3 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == 0 .and. & 
IEOR(s_a4,s_a2) == IEOR(s_a3,s_a1) .and. &
IEOR(IEOR(s_i,s_a0),s_a4) == IEOR(IEOR(s_a2,s_a3),s_a1)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D3(i,a0,a4,a2,a3,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a2, i_a3, i_a1) =  &
  D3_(s_a1, s_a3, s_a2, s_a4, s_a0, s_i)%array(i_a1, i_a3, i_a2, i_a4, i_a0, i_i)
end do
end do
end do
end do
end do
! Z2 <-- V2(a4,a2,a3,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a3, i_a1) =  &
  V2_(s_a1, s_a3, s_a2)%array(i_a1, i_a3, i_a2)
end do
end do
end do

! Z3 <-- W13(i,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W13(i,a0)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W13_(s_a0, s_i)%array(i_a0, i_i) = &
    W13_(s_a0, s_i)%array(i_a0, i_i) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no1_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no0_x0_type2_eri_o &
  (sb, ib, T2, W10, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: T2(*), W10(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xaa(sleft, W10, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no0_x0_type2_eri_o &
  (sb, ib, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no0_x0_type2_eri_o



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
subroutine g_sigma_covv_covv_no0_x0_type2_eri_o &
  (s_b, i_b, T2_, W10_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W10_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_a, i_a, s_i, i_i
! S2(w,i,a,b) += (    1.00000000) T2(w,a0,a,b) W10(i,a0) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_a,s_b) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(w,a0,a,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  T2_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w)
end do
end do
end do
! Z2 <-- W10(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  W10_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no0_x0_type2_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no1_x0_type2_eri_o &
  (sb, ib, T2, W13, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: T2(*), W13(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xaa(sleft, W13, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no1_x0_type2_eri_o &
  (sb, ib, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no1_x0_type2_eri_o



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
subroutine g_sigma_covv_covv_no1_x0_type2_eri_o &
  (s_b, i_b, T2_, W13_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W13_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_a, i_a, s_i, i_i
! S2(w,i,a,b) += (   -0.50000000) T2(a0,w,a,b) W13(i,a0) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a,s_b) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,w,a,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  T2_(s_a, s_w, s_a0)%array(i_a, i_w, i_a0)
end do
end do
end do
! Z2 <-- W13(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  W13_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 0.50000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no1_x0_type2_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no0_x0_type1_eri_v &
  (sv0, iv0, V2, W15, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: V2(*), W15(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xaav(sleft, W15, nir, nsym, psym) ! -> Xaav (allocate) 
call g_sigma_covv_covv_no0_x0_type1_eri_v &
  (sv0, iv0, h2_i, Xaav, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaav)

end subroutine g_if_sigma_covv_covv_no0_x0_type1_eri_v



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
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no0_x0_type1_eri_v &
  (s_v0, i_v0, V2_, W15_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W15_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a2, i_a2, s_a, i_a, s_i, i_i, s_a0, i_a0
! W15(i,a0,v0,a) += (    1.00000000) V2(v0,a1,a2,a) D2(i,a0,a1,a2) 
do s_a1 = 0, nir-1
do s_a2 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_v0,s_a) .and. & 
IEOR(s_v0,s_a1) == IEOR(s_a2,s_a) .and. &
IEOR(s_i,s_a0) == IEOR(s_a1,s_a2)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(v0,a1,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a1, i_a2) =  &
  V2_(s_a, s_a2, s_a1)%array(i_a, i_a2, i_a1)
end do
end do
end do
! Z2 <-- D2(i,a0,a1,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a2, i_i, i_a0) =  &
  D2_(s_a2, s_a1, s_a0, s_i)%array(i_a2, i_a1, i_a0, i_i)
end do
end do
end do
end do

! Z3 <-- W15(i,a0,v0,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W15(i,a0,v0,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W15_(s_a, s_a0, s_i)%array(i_a, i_a0, i_i) = &
    W15_(s_a, s_a0, s_i)%array(i_a, i_a0, i_i) &
  + Z3_(i_a, i_i, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no0_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no0_x1_type1_eri_v &
  (sb, ib, sv0, iv0, T2, W15, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W15(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xaav(sleft, W15, nir, nsym, psym) ! -> Xaav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no0_x1_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, Xaav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaav)

end subroutine g_if_sigma_covv_covv_no0_x1_type1_eri_v



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
subroutine g_sigma_covv_covv_no0_x1_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, W15_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W15_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_i, i_i, s_a, i_a
! S2(w,i,a,b) += (   -1.00000000) T2(w,a0,v0,b) W15(i,a0,v0,a) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_v0,s_b) .and. &
IEOR(s_i,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W15(i,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a, i_a0) =  &
  W15_(s_a, s_a0, s_i)%array(i_a, i_a0, i_i)
end do
end do
end do
! Z2 <-- T2(w,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  T2_(s_v0, s_a0, s_w)%array(i_v0, i_a0, i_w)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_i, i_a, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no0_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no1_x0_type1_eri_v &
  (sa, ia, sb, ib, T2, V2, W16, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sb, ib
real(kind=8), intent(inout) :: T2(*), V2(*), W16(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa,sb)

call set_symblock_Xca(sleft, W16, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_covv_covv_no1_x0_type1_eri_v &
  (sa, ia, sb, ib, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_covv_covv_no1_x0_type1_eri_v



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
subroutine g_sigma_covv_covv_no1_x0_type1_eri_v &
  (s_a, i_a, s_b, i_b, T2_, V2_, W16_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W16_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_v0, i_v0, s_c0, i_c0, s_a0, i_a0
! W16(w,a0,a,b) += (    1.00000000) V2(b,w,v0,c0) T2(c0,a0,v0,a) 
do s_w = 0, nir-1
do s_v0 = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_a,s_b) .and. & 
IEOR(s_b,s_w) == IEOR(s_v0,s_c0) .and. &
IEOR(s_c0,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- T2(c0,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_c0, i_v0) =  &
  T2_(s_v0, s_a0, s_c0)%array(i_v0, i_a0, i_c0)
end do
end do
end do
! Z2 <-- V2(b,w,v0,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_v0, i_w) =  &
  V2_(s_c0, s_v0, s_w)%array(i_c0, i_v0, i_w)
end do
end do
end do

! Z3 <-- W16(w,a0,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W16(w,a0,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W16_(s_a0, s_w)%array(i_a0, i_w) = &
    W16_(s_a0, s_w)%array(i_a0, i_w) &
  + Z3_(i_a0, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no1_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no1_x1_type1_eri_v &
  (sa, ia, sb, ib, W16, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sb, ib
real(kind=8), intent(inout) :: W16(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa,sb)

call set_symblock_Xca(sleft, W16, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no1_x1_type1_eri_v &
  (sa, ia, sb, ib, Xca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_covv_covv_no1_x1_type1_eri_v



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
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
subroutine g_sigma_covv_covv_no1_x1_type1_eri_v &
  (s_a, i_a, s_b, i_b, W16_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W16_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_w, i_w
! S2(w,i,a,b) += (   -2.00000000) D1(i,a0) W16(w,a0,a,b) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_a,s_b)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D1(i,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do
! Z2 <-- W16(w,a0,a,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  W16_(s_a0, s_w)%array(i_a0, i_w)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_i, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no1_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no2_x0_type1_eri_v &
  (sb, ib, sv0, iv0, V2, W17, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W17(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xaa(sleft, W17, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_covv_covv_no2_x0_type1_eri_v &
  (sb, ib, sv0, iv0, h2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no2_x0_type1_eri_v



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
subroutine g_sigma_covv_covv_no2_x0_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, V2_, W17_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W17_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W17(i,a0,v0,b) += (    1.00000000) V2(v0,a2,b,a1) D2(i,a1,a2,a0) 
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_v0,s_b) .and. & 
IEOR(s_v0,s_a2) == IEOR(s_b,s_a1) .and. &
IEOR(s_i,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- D2(i,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a1, i_a2) =  &
  D2_(s_a0, s_a2, s_a1, s_i)%array(i_a0, i_a2, i_a1, i_i)
end do
end do
end do
end do
! Z2 <-- V2(v0,a2,b,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a2) =  &
  V2_(s_a1, s_b, s_a2)%array(i_a1, i_b, i_a2)
end do
end do

! Z3 <-- W17(i,a0,v0,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W17(i,a0,v0,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W17_(s_a0, s_i)%array(i_a0, i_i) = &
    W17_(s_a0, s_i)%array(i_a0, i_i) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no2_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no2_x1_type1_eri_v &
  (sb, ib, sv0, iv0, T2, W17, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W17(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xaa(sleft, W17, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no2_x1_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no2_x1_type1_eri_v



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
subroutine g_sigma_covv_covv_no2_x1_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, W17_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W17_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_a, i_a, s_i, i_i
! S2(w,i,a,b) += (   -1.00000000) T2(a0,w,a,v0) W17(i,a0,v0,b) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a,s_v0) .and. &
IEOR(s_i,s_a0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,w,a,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  T2_(s_a, s_w, s_a0)%array(i_a, i_w, i_a0)
end do
end do
end do
! Z2 <-- W17(i,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  W17_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no2_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no3_x0_type1_eri_v &
  (sa, ia, sb, ib, T2, V2, W18, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sb, ib
real(kind=8), intent(inout) :: T2(*), V2(*), W18(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa,sb)

call set_symblock_Xca(sleft, W18, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_covv_covv_no3_x0_type1_eri_v &
  (sa, ia, sb, ib, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_covv_covv_no3_x0_type1_eri_v



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
subroutine g_sigma_covv_covv_no3_x0_type1_eri_v &
  (s_a, i_a, s_b, i_b, T2_, V2_, W18_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W18_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_w, i_w, s_c0, i_c0, s_a0, i_a0
! W18(w,a0,a,b) += (    1.00000000) V2(b,v0,w,c0) T2(c0,a0,v0,a) 
do s_v0 = 0, nir-1
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_a,s_b) .and. & 
IEOR(s_b,s_v0) == IEOR(s_w,s_c0) .and. &
IEOR(s_c0,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- T2(c0,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_c0, i_v0) =  &
  T2_(s_v0, s_a0, s_c0)%array(i_v0, i_a0, i_c0)
end do
end do
end do
! Z2 <-- V2(b,v0,w,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_v0, i_w) =  &
  V2_(s_c0, s_w, s_v0)%array(i_c0, i_w, i_v0)
end do
end do
end do

! Z3 <-- W18(w,a0,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W18(w,a0,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W18_(s_a0, s_w)%array(i_a0, i_w) = &
    W18_(s_a0, s_w)%array(i_a0, i_w) &
  + Z3_(i_a0, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no3_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no3_x1_type1_eri_v &
  (sa, ia, sb, ib, W18, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sb, ib
real(kind=8), intent(inout) :: W18(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa,sb)

call set_symblock_Xca(sleft, W18, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no3_x1_type1_eri_v &
  (sa, ia, sb, ib, Xca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_covv_covv_no3_x1_type1_eri_v



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
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
subroutine g_sigma_covv_covv_no3_x1_type1_eri_v &
  (s_a, i_a, s_b, i_b, W18_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W18_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_w, i_w
! S2(w,i,a,b) += (    1.00000000) D1(i,a0) W18(w,a0,a,b) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_a,s_b)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D1(i,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do
! Z2 <-- W18(w,a0,a,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  W18_(s_a0, s_w)%array(i_a0, i_w)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_i, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no3_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no4_x0_type1_eri_v &
  (sb, ib, sv0, iv0, V2, W19, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W19(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xaa(sleft, W19, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_covv_covv_no4_x0_type1_eri_v &
  (sb, ib, sv0, iv0, h2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no4_x0_type1_eri_v



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
subroutine g_sigma_covv_covv_no4_x0_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, V2_, W19_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W19_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W19(i,a0,v0,b) += (    1.00000000) V2(v0,b,a2,a1) D2(i,a0,a2,a1) 
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_v0,s_b) .and. & 
IEOR(s_v0,s_b) == IEOR(s_a2,s_a1) .and. &
IEOR(s_i,s_a0) == IEOR(s_a2,s_a1)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(i,a0,a2,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a2, i_a1) =  &
  D2_(s_a1, s_a2, s_a0, s_i)%array(i_a1, i_a2, i_a0, i_i)
end do
end do
end do
end do
! Z2 <-- V2(v0,b,a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1) =  &
  V2_(s_a1, s_a2, s_b)%array(i_a1, i_a2, i_b)
end do
end do

! Z3 <-- W19(i,a0,v0,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W19(i,a0,v0,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W19_(s_a0, s_i)%array(i_a0, i_i) = &
    W19_(s_a0, s_i)%array(i_a0, i_i) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no4_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no4_x1_type1_eri_v &
  (sb, ib, sv0, iv0, T2, W19, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W19(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xaa(sleft, W19, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no4_x1_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no4_x1_type1_eri_v



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
subroutine g_sigma_covv_covv_no4_x1_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, W19_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W19_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_a, i_a, s_i, i_i
! S2(w,i,a,b) += (   -1.00000000) T2(a0,w,a,v0) W19(i,a0,v0,b) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a,s_v0) .and. &
IEOR(s_i,s_a0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,w,a,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  T2_(s_a, s_w, s_a0)%array(i_a, i_w, i_a0)
end do
end do
end do
! Z2 <-- W19(i,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  W19_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no4_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no5_x0_type1_eri_v &
  (sv0, iv0, V2, W21, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: V2(*), W21(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xaav(sleft, W21, nir, nsym, psym) ! -> Xaav (allocate) 
call g_sigma_covv_covv_no5_x0_type1_eri_v &
  (sv0, iv0, h2_i, Xaav, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaav)

end subroutine g_if_sigma_covv_covv_no5_x0_type1_eri_v



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
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no5_x0_type1_eri_v &
  (s_v0, i_v0, V2_, W21_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W21_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_a2, i_a2, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W21(i,a0,v0,a) += (    1.00000000) V2(v0,a,a2,a1) D2(i,a0,a2,a1) 
do s_a = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_v0,s_a) .and. & 
IEOR(s_v0,s_a) == IEOR(s_a2,s_a1) .and. &
IEOR(s_i,s_a0) == IEOR(s_a2,s_a1)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(v0,a,a2,a1) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a2, i_a1) =  &
  V2_(s_a1, s_a2, s_a)%array(i_a1, i_a2, i_a)
end do
end do
end do
! Z2 <-- D2(i,a0,a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1, i_i, i_a0) =  &
  D2_(s_a1, s_a2, s_a0, s_i)%array(i_a1, i_a2, i_a0, i_i)
end do
end do
end do
end do

! Z3 <-- W21(i,a0,v0,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W21(i,a0,v0,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W21_(s_a, s_a0, s_i)%array(i_a, i_a0, i_i) = &
    W21_(s_a, s_a0, s_i)%array(i_a, i_a0, i_i) &
  + Z3_(i_a, i_i, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no5_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no5_x1_type1_eri_v &
  (sb, ib, sv0, iv0, T2, W21, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W21(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xaav(sleft, W21, nir, nsym, psym) ! -> Xaav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no5_x1_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, Xaav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaav)

end subroutine g_if_sigma_covv_covv_no5_x1_type1_eri_v



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
subroutine g_sigma_covv_covv_no5_x1_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, W21_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W21_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_i, i_i, s_a, i_a
! S2(w,i,a,b) += (    2.00000000) T2(w,a0,v0,b) W21(i,a0,v0,a) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_v0,s_b) .and. &
IEOR(s_i,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W21(i,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a, i_a0) =  &
  W21_(s_a, s_a0, s_i)%array(i_a, i_a0, i_i)
end do
end do
end do
! Z2 <-- T2(w,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  T2_(s_v0, s_a0, s_w)%array(i_v0, i_a0, i_w)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_i, i_a, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no5_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no6_x0_type1_eri_v &
  (sv0, iv0, V2, W23, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: V2(*), W23(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xaav(sleft, W23, nir, nsym, psym) ! -> Xaav (allocate) 
call g_sigma_covv_covv_no6_x0_type1_eri_v &
  (sv0, iv0, h2_i, Xaav, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaav)

end subroutine g_if_sigma_covv_covv_no6_x0_type1_eri_v



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
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no6_x0_type1_eri_v &
  (s_v0, i_v0, V2_, W23_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W23_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a2, i_a2, s_a, i_a, s_i, i_i, s_a0, i_a0
! W23(i,a0,v0,a) += (    1.00000000) V2(v0,a1,a2,a) D2(i,a2,a1,a0) 
do s_a1 = 0, nir-1
do s_a2 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_v0,s_a) .and. & 
IEOR(s_v0,s_a1) == IEOR(s_a2,s_a) .and. &
IEOR(s_i,s_a2) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(v0,a1,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a1, i_a2) =  &
  V2_(s_a, s_a2, s_a1)%array(i_a, i_a2, i_a1)
end do
end do
end do
! Z2 <-- D2(i,a2,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a2, i_i, i_a0) =  &
  D2_(s_a0, s_a1, s_a2, s_i)%array(i_a0, i_a1, i_a2, i_i)
end do
end do
end do
end do

! Z3 <-- W23(i,a0,v0,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W23(i,a0,v0,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W23_(s_a, s_a0, s_i)%array(i_a, i_a0, i_i) = &
    W23_(s_a, s_a0, s_i)%array(i_a, i_a0, i_i) &
  + Z3_(i_a, i_i, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no6_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no6_x1_type1_eri_v &
  (sb, ib, sv0, iv0, T2, W23, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W23(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xaav(sleft, W23, nir, nsym, psym) ! -> Xaav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no6_x1_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, Xaav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaav)

end subroutine g_if_sigma_covv_covv_no6_x1_type1_eri_v



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
subroutine g_sigma_covv_covv_no6_x1_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, W23_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W23_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_i, i_i, s_a, i_a
! S2(w,i,a,b) += (   -1.00000000) T2(w,a0,b,v0) W23(i,a0,v0,a) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_b,s_v0) .and. &
IEOR(s_i,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W23(i,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a, i_a0) =  &
  W23_(s_a, s_a0, s_i)%array(i_a, i_a0, i_i)
end do
end do
end do
! Z2 <-- T2(w,a0,b,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  T2_(s_b, s_a0, s_w)%array(i_b, i_a0, i_w)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_i, i_a, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no6_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no7_x0_type1_eri_v &
  (sb, ib, sv0, iv0, T2, V2, W24, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W24(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xcav(sleft, W24, nir, nsym, psym) ! -> Xcav (allocate) 
call g_sigma_covv_covv_no7_x0_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, h2_i, Xcav, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no7_x0_type1_eri_v



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
subroutine g_sigma_covv_covv_no7_x0_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, V2_, W24_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W24_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_a0, i_a0, s_a, i_a
! W24(w,a0,a,b) += (    1.00000000) V2(b,w,v0,c0) T2(c0,a0,a,v0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_a,s_b) .and. & 
IEOR(s_b,s_w) == IEOR(s_v0,s_c0) .and. &
IEOR(s_c0,s_a0) == IEOR(s_a,s_v0)) then

if(psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c0,a0,a,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a, i_c0) =  &
  T2_(s_a, s_a0, s_c0)%array(i_a, i_a0, i_c0)
end do
end do
end do
! Z2 <-- V2(b,w,v0,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  V2_(s_c0, s_v0, s_w)%array(i_c0, i_v0, i_w)
end do
end do

! Z3 <-- W24(w,a0,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_a))

! W24(w,a0,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W24_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w) = &
    W24_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w) &
  + Z3_(i_a0, i_a, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no7_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no7_x1_type1_eri_v &
  (sb, ib, W24, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: W24(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sb

call set_symblock_Xcav(sleft, W24, nir, nsym, psym) ! -> Xcav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no7_x1_type1_eri_v &
  (sb, ib, Xcav, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no7_x1_type1_eri_v



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
!    -- is4RDM.second true
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no7_x1_type1_eri_v &
  (s_b, i_b, W24_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W24_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_w, i_w, s_a, i_a
! S2(w,i,a,b) += (    1.00000000) D1(i,a0) W24(w,a0,a,b) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_a,s_b)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W24(w,a0,a,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  W24_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w)
end do
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no7_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no8_x0_type1_eri_v &
  (sb, ib, sv0, iv0, V2, W25, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W25(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xaa(sleft, W25, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_covv_covv_no8_x0_type1_eri_v &
  (sb, ib, sv0, iv0, h2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no8_x0_type1_eri_v



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
subroutine g_sigma_covv_covv_no8_x0_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, V2_, W25_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W25_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W25(i,a0,v0,b) += (    1.00000000) V2(v0,a2,b,a1) D2(i,a1,a2,a0) 
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_v0,s_b) .and. & 
IEOR(s_v0,s_a2) == IEOR(s_b,s_a1) .and. &
IEOR(s_i,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- D2(i,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a1, i_a2) =  &
  D2_(s_a0, s_a2, s_a1, s_i)%array(i_a0, i_a2, i_a1, i_i)
end do
end do
end do
end do
! Z2 <-- V2(v0,a2,b,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a2) =  &
  V2_(s_a1, s_b, s_a2)%array(i_a1, i_b, i_a2)
end do
end do

! Z3 <-- W25(i,a0,v0,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W25(i,a0,v0,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W25_(s_a0, s_i)%array(i_a0, i_i) = &
    W25_(s_a0, s_i)%array(i_a0, i_i) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no8_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no8_x1_type1_eri_v &
  (sb, ib, sv0, iv0, T2, W25, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W25(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xaa(sleft, W25, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no8_x1_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no8_x1_type1_eri_v



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
subroutine g_sigma_covv_covv_no8_x1_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, W25_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W25_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_a, i_a, s_i, i_i
! S2(w,i,a,b) += (    2.00000000) T2(w,a0,a,v0) W25(i,a0,v0,b) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_a,s_v0) .and. &
IEOR(s_i,s_a0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(w,a0,a,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  T2_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w)
end do
end do
end do
! Z2 <-- W25(i,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  W25_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no8_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no9_x0_type1_eri_v &
  (sb, ib, sv0, iv0, T2, V2, W26, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W26(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xcav(sleft, W26, nir, nsym, psym) ! -> Xcav (allocate) 
call g_sigma_covv_covv_no9_x0_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, h2_i, Xcav, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no9_x0_type1_eri_v



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
subroutine g_sigma_covv_covv_no9_x0_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, V2_, W26_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W26_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_a0, i_a0, s_a, i_a
! W26(w,a0,a,b) += (    1.00000000) V2(b,v0,w,c0) T2(c0,a0,a,v0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_a,s_b) .and. & 
IEOR(s_b,s_v0) == IEOR(s_w,s_c0) .and. &
IEOR(s_c0,s_a0) == IEOR(s_a,s_v0)) then

if(psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c0,a0,a,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a, i_c0) =  &
  T2_(s_a, s_a0, s_c0)%array(i_a, i_a0, i_c0)
end do
end do
end do
! Z2 <-- V2(b,v0,w,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  V2_(s_c0, s_w, s_v0)%array(i_c0, i_w, i_v0)
end do
end do

! Z3 <-- W26(w,a0,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_a))

! W26(w,a0,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W26_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w) = &
    W26_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w) &
  + Z3_(i_a0, i_a, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no9_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no9_x1_type1_eri_v &
  (sb, ib, W26, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: W26(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sb

call set_symblock_Xcav(sleft, W26, nir, nsym, psym) ! -> Xcav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no9_x1_type1_eri_v &
  (sb, ib, Xcav, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no9_x1_type1_eri_v



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
!    -- is4RDM.second true
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no9_x1_type1_eri_v &
  (s_b, i_b, W26_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W26_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_w, i_w, s_a, i_a
! S2(w,i,a,b) += (   -2.00000000) D1(i,a0) W26(w,a0,a,b) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_a,s_b)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W26(w,a0,a,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  W26_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w)
end do
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no9_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no10_x0_type1_eri_v &
  (sb, ib, sv0, iv0, V2, W27, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W27(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xaa(sleft, W27, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_covv_covv_no10_x0_type1_eri_v &
  (sb, ib, sv0, iv0, h2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no10_x0_type1_eri_v



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
subroutine g_sigma_covv_covv_no10_x0_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, V2_, W27_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W27_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W27(i,a0,v0,b) += (    1.00000000) V2(v0,b,a2,a1) D2(i,a0,a2,a1) 
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_v0,s_b) .and. & 
IEOR(s_v0,s_b) == IEOR(s_a2,s_a1) .and. &
IEOR(s_i,s_a0) == IEOR(s_a2,s_a1)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(i,a0,a2,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a2, i_a1) =  &
  D2_(s_a1, s_a2, s_a0, s_i)%array(i_a1, i_a2, i_a0, i_i)
end do
end do
end do
end do
! Z2 <-- V2(v0,b,a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1) =  &
  V2_(s_a1, s_a2, s_b)%array(i_a1, i_a2, i_b)
end do
end do

! Z3 <-- W27(i,a0,v0,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W27(i,a0,v0,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W27_(s_a0, s_i)%array(i_a0, i_i) = &
    W27_(s_a0, s_i)%array(i_a0, i_i) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no10_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no10_x1_type1_eri_v &
  (sb, ib, sv0, iv0, T2, W27, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W27(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xaa(sleft, W27, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no10_x1_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_covv_covv_no10_x1_type1_eri_v



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
subroutine g_sigma_covv_covv_no10_x1_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, W27_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W27_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_a, i_a, s_i, i_i
! S2(w,i,a,b) += (    2.00000000) T2(w,a0,a,v0) W27(i,a0,v0,b) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_a,s_v0) .and. &
IEOR(s_i,s_a0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(w,a0,a,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  T2_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w)
end do
end do
end do
! Z2 <-- W27(i,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  W27_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no10_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no11_x0_type1_eri_v &
  (sv0, iv0, V2, W29, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: V2(*), W29(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xaav(sleft, W29, nir, nsym, psym) ! -> Xaav (allocate) 
call g_sigma_covv_covv_no11_x0_type1_eri_v &
  (sv0, iv0, h2_i, Xaav, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaav)

end subroutine g_if_sigma_covv_covv_no11_x0_type1_eri_v



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
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no11_x0_type1_eri_v &
  (s_v0, i_v0, V2_, W29_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W29_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_a2, i_a2, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W29(i,a0,v0,a) += (    1.00000000) V2(v0,a,a2,a1) D2(i,a0,a2,a1) 
do s_a = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_v0,s_a) .and. & 
IEOR(s_v0,s_a) == IEOR(s_a2,s_a1) .and. &
IEOR(s_i,s_a0) == IEOR(s_a2,s_a1)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(v0,a,a2,a1) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a2, i_a1) =  &
  V2_(s_a1, s_a2, s_a)%array(i_a1, i_a2, i_a)
end do
end do
end do
! Z2 <-- D2(i,a0,a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1, i_i, i_a0) =  &
  D2_(s_a1, s_a2, s_a0, s_i)%array(i_a1, i_a2, i_a0, i_i)
end do
end do
end do
end do

! Z3 <-- W29(i,a0,v0,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W29(i,a0,v0,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W29_(s_a, s_a0, s_i)%array(i_a, i_a0, i_i) = &
    W29_(s_a, s_a0, s_i)%array(i_a, i_a0, i_i) &
  + Z3_(i_a, i_i, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no11_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no11_x1_type1_eri_v &
  (sb, ib, sv0, iv0, T2, W29, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W29(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xaav(sleft, W29, nir, nsym, psym) ! -> Xaav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no11_x1_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, Xaav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaav)

end subroutine g_if_sigma_covv_covv_no11_x1_type1_eri_v



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
subroutine g_sigma_covv_covv_no11_x1_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, W29_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W29_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_i, i_i, s_a, i_a
! S2(w,i,a,b) += (   -1.00000000) T2(w,a0,b,v0) W29(i,a0,v0,a) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_b,s_v0) .and. &
IEOR(s_i,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W29(i,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a, i_a0) =  &
  W29_(s_a, s_a0, s_i)%array(i_a, i_a0, i_i)
end do
end do
end do
! Z2 <-- T2(w,a0,b,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  T2_(s_b, s_a0, s_w)%array(i_b, i_a0, i_w)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_i, i_a, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no11_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no12_x0_type1_eri_v &
  (sb, ib, sv1, iv1, T2, V2, W30, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), W30(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xcav(sleft, W30, nir, nsym, psym) ! -> Xcav (allocate) 
call g_sigma_covv_covv_no12_x0_type1_eri_v &
  (sb, ib, sv1, iv1, av2_i, h2_i, Xcav, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no12_x0_type1_eri_v



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
subroutine g_sigma_covv_covv_no12_x0_type1_eri_v &
  (s_b, i_b, s_v1, i_v1, T2_, V2_, W30_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W30_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_a, i_a, s_w, i_w, s_a0, i_a0
! W30(w,a0,b,a) += (    1.00000000) V2(b,v1,v0,a) T2(w,a0,v0,v1) 
do s_v0 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_b,s_a) .and. & 
IEOR(s_b,s_v1) == IEOR(s_v0,s_a) .and. &
IEOR(s_w,s_a0) == IEOR(s_v0,s_v1)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- V2(b,v1,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_v0) =  &
  V2_(s_a, s_v0, s_v1)%array(i_a, i_v0, i_v1)
end do
end do
! Z2 <-- T2(w,a0,v0,v1) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_w, i_a0) =  &
  T2_(s_v0, s_a0, s_w)%array(i_v0, i_a0, i_w)
end do
end do
end do

! Z3 <-- W30(w,a0,b,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W30(w,a0,b,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W30_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w) = &
    W30_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w) &
  + Z3_(i_a, i_w, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no12_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no12_x1_type1_eri_v &
  (sb, ib, W30, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: W30(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sb

call set_symblock_Xcav(sleft, W30, nir, nsym, psym) ! -> Xcav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no12_x1_type1_eri_v &
  (sb, ib, Xcav, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no12_x1_type1_eri_v



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
!    -- is4RDM.second true
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no12_x1_type1_eri_v &
  (s_b, i_b, W30_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W30_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_w, i_w, s_a, i_a
! S2(w,i,a,b) += (    2.00000000) D1(i,a0) W30(w,a0,b,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_b,s_a)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W30(w,a0,b,a) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  W30_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w)
end do
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no12_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no13_x0_type1_eri_v &
  (sb, ib, sv0, iv0, T2, V2, W31, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W31(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xcav(sleft, W31, nir, nsym, psym) ! -> Xcav (allocate) 
call g_sigma_covv_covv_no13_x0_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, h2_i, Xcav, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no13_x0_type1_eri_v



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
subroutine g_sigma_covv_covv_no13_x0_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, V2_, W31_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W31_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v1, i_v1, s_a, i_a, s_w, i_w, s_a0, i_a0
! W31(w,a0,b,a) += (    1.00000000) V2(b,v1,v0,a) T2(w,a0,v1,v0) 
do s_v1 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_b,s_a) .and. & 
IEOR(s_b,s_v1) == IEOR(s_v0,s_a) .and. &
IEOR(s_w,s_a0) == IEOR(s_v1,s_v0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_V, s_v1) > 0) then

! Z1 <-- V2(b,v1,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v1):psym(I_END,I_V, s_v1)))
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_v1) =  &
  V2_(s_a, s_v0, s_v1)%array(i_a, i_v0, i_v1)
end do
end do
! Z2 <-- T2(w,a0,v1,v0) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v1):psym(I_END,I_V, s_v1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Z2_(i_v1, i_w, i_a0) =  &
  T2_(s_v1, s_a0, s_w)%array(i_v1, i_a0, i_w)
end do
end do
end do

! Z3 <-- W31(w,a0,b,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_V, s_v1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W31(w,a0,b,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W31_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w) = &
    W31_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w) &
  + Z3_(i_a, i_w, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_V, s_v1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no13_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_covv_covv_no13_x1_type1_eri_v &
  (sb, ib, W31, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: W31(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sb

call set_symblock_Xcav(sleft, W31, nir, nsym, psym) ! -> Xcav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_covv_covv_no13_x1_type1_eri_v &
  (sb, ib, Xcav, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcav)

end subroutine g_if_sigma_covv_covv_no13_x1_type1_eri_v



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
!    -- is4RDM.second true
!    -- allRDM.second true
!    -- precedence    -1
subroutine g_sigma_covv_covv_no13_x1_type1_eri_v &
  (s_b, i_b, W31_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W31_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_w, i_w, s_a, i_a
! S2(w,i,a,b) += (   -1.00000000) D1(i,a0) W31(w,a0,b,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a,s_b) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_b,s_a)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W31(w,a0,b,a) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a, i_a0) =  &
  W31_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w)
end do
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,i,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) = &
    S2_(s_a, s_i, s_w)%array(i_a, i_i, i_w) &
  + Z3_(i_w, i_a, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_covv_covv_no13_x1_type1_eri_v

