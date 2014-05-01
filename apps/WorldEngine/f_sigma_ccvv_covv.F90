#include <sci/icmr/fsrc/f_mr.fh>


!  `7MM"""YMM                         mm               
!    MM    `7                         MM                  
!    MM   d  .gP"Ya `7MMpMMMb.pMMMb.mmMMmm ,pW"Wq.      
!    MM""MM ,M'   Yb  MM    MM    MM  MM  6W'   `Wb   
!    MM   Y 8M""""""  MM    MM    MM  MM  8M     M8 
!    MM     YM.    ,  MM    MM    MM  MM  YA.   ,A9       
!  .JMML.    `Mbmmd'.JMML  JMML  JMML.`Mbmo`Ybmd9'        

!                                    Generated date : Sun Apr 20 10:26:25 2014



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no0_x0_type0_noeri &
  (W0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

real(kind=8), intent(inout) :: W0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W0, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_ccvv_covv_no0_x0_type0_noeri &
  (Xca, d1, fc1, nir, nsym, psym, flops)

deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no0_x0_type0_noeri



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
subroutine g_sigma_ccvv_covv_no0_x0_type0_noeri &
  (W0_, D1_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W0_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a0, i_a0, s_w, i_w
! W0(w,a0) += (    1.00000000) D1(a1,a0) Fc1(w,a1) 
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_a0) == 0 .and. & 
IEOR(s_a1,s_a0) == 0 .and. &
IEOR(s_w,s_a1) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do
! Z2 <-- Fc1(w,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w) =  &
  Fc1_(s_a1, s_w)%array(i_a1, i_w)
end do
end do

! Z3 <-- W0(w,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W0(w,a0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W0_(s_a0, s_w)%array(i_a0, i_w) = &
    W0_(s_a0, s_w)%array(i_a0, i_w) &
  + Z3_(i_a0, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no0_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no0_x1_type0_noeri &
  (sb, ib, T2, W0, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: T2(*), W0(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W0, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no0_x1_type0_noeri &
  (sb, ib, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no0_x1_type0_noeri



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
subroutine g_sigma_ccvv_covv_no0_x1_type0_noeri &
  (s_b, i_b, T2_, W0_, S2_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W0_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_a0, i_a0, s_a, i_a, s_w, i_w
! S2(w,x,a,b) += (    1.00000000) T2(x,a0,a,b) W0(w,a0) 
do s_x = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_x,s_a0) == IEOR(s_a,s_b) .and. &
IEOR(s_w,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(x,a0,a,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_a, i_a0) =  &
  T2_(s_a, s_a0, s_x)%array(i_a, i_a0, i_x)
end do
end do
end do
! Z2 <-- W0(w,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  W0_(s_a0, s_w)%array(i_a0, i_w)
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_x, i_a, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a) * &
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

end subroutine g_sigma_ccvv_covv_no0_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no1_x0_type0_noeri &
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
call set_symblock_Xca(sleft, W1, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_ccvv_covv_no1_x0_type0_noeri &
  (Xca, d1, fc1, nir, nsym, psym, flops)

deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no1_x0_type0_noeri



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
subroutine g_sigma_ccvv_covv_no1_x0_type0_noeri &
  (W1_, D1_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W1_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a0, i_a0, s_x, i_x
! W1(x,a0) += (    1.00000000) D1(a1,a0) Fc1(x,a1) 
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_x,s_a0) == 0 .and. & 
IEOR(s_a1,s_a0) == 0 .and. &
IEOR(s_x,s_a1) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do
! Z2 <-- Fc1(x,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_x) =  &
  Fc1_(s_a1, s_x)%array(i_a1, i_x)
end do
end do

! Z3 <-- W1(x,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W1(x,a0)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W1_(s_a0, s_x)%array(i_a0, i_x) = &
    W1_(s_a0, s_x)%array(i_a0, i_x) &
  + Z3_(i_a0, i_x)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no1_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no1_x1_type0_noeri &
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
call set_symblock_Xca(sleft, W1, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no1_x1_type0_noeri &
  (sb, ib, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no1_x1_type0_noeri



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
subroutine g_sigma_ccvv_covv_no1_x1_type0_noeri &
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

integer :: s_w, i_w, s_a0, i_a0, s_a, i_a, s_x, i_x
! S2(w,x,a,b) += (   -2.00000000) T2(w,a0,a,b) W1(x,a0) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_a,s_b) .and. &
IEOR(s_x,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
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
! Z2 <-- W1(x,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x) =  &
  W1_(s_a0, s_x)%array(i_a0, i_x)
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_w, i_a, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no1_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no2_x0_type0_noeri &
  (W2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

real(kind=8), intent(inout) :: W2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W2, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_ccvv_covv_no2_x0_type0_noeri &
  (Xca, d1, fc1, nir, nsym, psym, flops)

deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no2_x0_type0_noeri



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
subroutine g_sigma_ccvv_covv_no2_x0_type0_noeri &
  (W2_, D1_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W2_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a0, i_a0, s_w, i_w
! W2(w,a0) += (    1.00000000) D1(a1,a0) Fc1(w,a1) 
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_a0) == 0 .and. & 
IEOR(s_a1,s_a0) == 0 .and. &
IEOR(s_w,s_a1) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do
! Z2 <-- Fc1(w,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w) =  &
  Fc1_(s_a1, s_w)%array(i_a1, i_w)
end do
end do

! Z3 <-- W2(w,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W2(w,a0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W2_(s_a0, s_w)%array(i_a0, i_w) = &
    W2_(s_a0, s_w)%array(i_a0, i_w) &
  + Z3_(i_a0, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no2_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no2_x1_type0_noeri &
  (sb, ib, T2, W2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: T2(*), W2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W2, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no2_x1_type0_noeri &
  (sb, ib, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no2_x1_type0_noeri



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
subroutine g_sigma_ccvv_covv_no2_x1_type0_noeri &
  (s_b, i_b, T2_, W2_, S2_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W2_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_x, i_x, s_a, i_a, s_w, i_w
! S2(w,x,a,b) += (   -2.00000000) T2(a0,x,a,b) W2(w,a0) 
do s_a0 = 0, nir-1
do s_x = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_x) == IEOR(s_a,s_b) .and. &
IEOR(s_w,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,x,a,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_a, i_a0) =  &
  T2_(s_a, s_x, s_a0)%array(i_a, i_x, i_a0)
end do
end do
end do
! Z2 <-- W2(w,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  W2_(s_a0, s_w)%array(i_a0, i_w)
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_x, i_a, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a) * &
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

end subroutine g_sigma_ccvv_covv_no2_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no3_x0_type0_noeri &
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
call set_symblock_Xca(sleft, W3, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_ccvv_covv_no3_x0_type0_noeri &
  (Xca, d1, fc1, nir, nsym, psym, flops)

deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no3_x0_type0_noeri



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
subroutine g_sigma_ccvv_covv_no3_x0_type0_noeri &
  (W3_, D1_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W3_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a0, i_a0, s_x, i_x
! W3(x,a0) += (    1.00000000) D1(a1,a0) Fc1(x,a1) 
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_x,s_a0) == 0 .and. & 
IEOR(s_a1,s_a0) == 0 .and. &
IEOR(s_x,s_a1) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do
! Z2 <-- Fc1(x,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_x) =  &
  Fc1_(s_a1, s_x)%array(i_a1, i_x)
end do
end do

! Z3 <-- W3(x,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W3(x,a0)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W3_(s_a0, s_x)%array(i_a0, i_x) = &
    W3_(s_a0, s_x)%array(i_a0, i_x) &
  + Z3_(i_a0, i_x)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no3_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no3_x1_type0_noeri &
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
call set_symblock_Xca(sleft, W3, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no3_x1_type0_noeri &
  (sb, ib, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no3_x1_type0_noeri



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
subroutine g_sigma_ccvv_covv_no3_x1_type0_noeri &
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

integer :: s_a0, i_a0, s_w, i_w, s_a, i_a, s_x, i_x
! S2(w,x,a,b) += (    1.00000000) T2(a0,w,a,b) W3(x,a0) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a,s_b) .and. &
IEOR(s_x,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
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
! Z2 <-- W3(x,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x) =  &
  W3_(s_a0, s_x)%array(i_a0, i_x)
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_w, i_a, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no3_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no0_x0_type1_eri_c &
  (sx, ix, V2, W4, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sx, ix
real(kind=8), intent(inout) :: V2(*), W4(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sx, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sx

call set_symblock_Xcca(sleft, W4, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_ccvv_covv_no0_x0_type1_eri_c &
  (sx, ix, h2_i, Xcca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccvv_covv_no0_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no0_x0_type1_eri_c &
  (s_x, i_x, V2_, W4_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W4_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_w, i_w, s_a1, i_a1, s_a0, i_a0
! W4(x,w,c0,a0) += (    1.00000000) V2(x,c0,w,a1) D1(a1,a0) 
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_w) == IEOR(s_c0,s_a0) .and. & 
IEOR(s_x,s_c0) == IEOR(s_w,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do
! Z2 <-- V2(x,c0,w,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_c0, i_w) =  &
  V2_(s_a1, s_w, s_c0)%array(i_a1, i_w, i_c0)
end do
end do
end do

! Z3 <-- W4(x,w,c0,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W4(x,w,c0,a0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W4_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w) = &
    W4_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w) &
  + Z3_(i_a0, i_c0, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no0_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no0_x1_type1_eri_c &
  (sb, ib, sx, ix, T2, W4, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sx, ix
real(kind=8), intent(inout) :: T2(*), W4(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sx

call set_symblock_Xcca(sleft, W4, nir, nsym, psym) ! -> Xcca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no0_x1_type1_eri_c &
  (sb, ib, sx, ix, av2_i, Xcca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccvv_covv_no0_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no0_x1_type1_eri_c &
  (s_b, i_b, s_x, i_x, T2_, W4_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W4_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a0, i_a0, s_a, i_a, s_w, i_w
! S2(w,x,a,b) += (   -1.00000000) T2(c0,a0,a,b) W4(x,w,c0,a0) 
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_c0,s_a0) == IEOR(s_a,s_b) .and. &
IEOR(s_x,s_w) == IEOR(s_c0,s_a0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
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
! Z2 <-- W4(x,w,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0, i_w) =  &
  W4_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_a, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no0_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no1_x0_type1_eri_c &
  (sx, ix, V2, W5, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sx, ix
real(kind=8), intent(inout) :: V2(*), W5(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sx, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sx

call set_symblock_Xcca(sleft, W5, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_ccvv_covv_no1_x0_type1_eri_c &
  (sx, ix, h2_i, Xcca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccvv_covv_no1_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no1_x0_type1_eri_c &
  (s_x, i_x, V2_, W5_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W5_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_w, i_w, s_c0, i_c0, s_a0, i_a0
! W5(x,w,c0,a0) += (    1.00000000) V2(x,a1,w,c0) D1(a1,a0) 
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_w) == IEOR(s_c0,s_a0) .and. & 
IEOR(s_x,s_a1) == IEOR(s_w,s_c0) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do
! Z2 <-- V2(x,a1,w,c0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w, i_c0) =  &
  V2_(s_c0, s_w, s_a1)%array(i_c0, i_w, i_a1)
end do
end do
end do

! Z3 <-- W5(x,w,c0,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W5(x,w,c0,a0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W5_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w) = &
    W5_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w) &
  + Z3_(i_a0, i_w, i_c0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no1_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no1_x1_type1_eri_c &
  (sb, ib, sx, ix, T2, W5, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sx, ix
real(kind=8), intent(inout) :: T2(*), W5(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sx

call set_symblock_Xcca(sleft, W5, nir, nsym, psym) ! -> Xcca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no1_x1_type1_eri_c &
  (sb, ib, sx, ix, av2_i, Xcca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccvv_covv_no1_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no1_x1_type1_eri_c &
  (s_b, i_b, s_x, i_x, T2_, W5_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W5_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a0, i_a0, s_a, i_a, s_w, i_w
! S2(w,x,a,b) += (    2.00000000) T2(c0,a0,a,b) W5(x,w,c0,a0) 
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_c0,s_a0) == IEOR(s_a,s_b) .and. &
IEOR(s_x,s_w) == IEOR(s_c0,s_a0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
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
! Z2 <-- W5(x,w,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0, i_w) =  &
  W5_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_a, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no1_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no2_x0_type1_eri_c &
  (sw, iw, V2, W6, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), W6(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_Xa(sleft, W6, nir, nsym, psym) ! -> Xa (allocate) 
call g_sigma_ccvv_covv_no2_x0_type1_eri_c &
  (sw, iw, h2_i, Xa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xa)

end subroutine g_if_sigma_ccvv_covv_no2_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no2_x0_type1_eri_c &
  (s_w, i_w, V2_, W6_, D2_, nir, nsym, psym, flops)

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
type(symblock1), intent(inout) :: W6_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a3, i_a3, s_a1, i_a1, s_a0, i_a0
! W6(w,a0) += (    1.00000000) V2(w,a2,a3,a1) D2(a3,a1,a2,a0) 
do s_a2 = 0, nir-1
do s_a3 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == 0 .and. & 
IEOR(s_w,s_a2) == IEOR(s_a3,s_a1) .and. &
IEOR(s_a3,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- D2(a3,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a3, i_a1, i_a2) =  &
  D2_(s_a0, s_a2, s_a1, s_a3)%array(i_a0, i_a2, i_a1, i_a3)
end do
end do
end do
end do
! Z2 <-- V2(w,a2,a3,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a1, i_a2) =  &
  V2_(s_a1, s_a3, s_a2)%array(i_a1, i_a3, i_a2)
end do
end do
end do

! Z3 <-- W6(w,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W6(w,a0)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W6_(s_a0)%array(i_a0) = &
    W6_(s_a0)%array(i_a0) &
  + Z3_(i_a0)
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no2_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no2_x1_type1_eri_c &
  (sb, ib, sw, iw, T2, W6, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: T2(*), W6(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xa(sleft, W6, nir, nsym, psym) ! -> Xa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no2_x1_type1_eri_c &
  (sb, ib, sw, iw, av2_i, Xa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xa)

end subroutine g_if_sigma_ccvv_covv_no2_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no2_x1_type1_eri_c &
  (s_b, i_b, s_w, i_w, T2_, W6_, S2_, nir, nsym, psym, flops)

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
type(symblock1), intent(inout) :: W6_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_a0, i_a0, s_a, i_a
! S2(w,x,a,b) += (    1.00000000) T2(x,a0,a,b) W6(w,a0) 
do s_x = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_x,s_a0) == IEOR(s_a,s_b) .and. &
IEOR(s_w,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(x,a0,a,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_a, i_a0) =  &
  T2_(s_a, s_a0, s_x)%array(i_a, i_a0, i_x)
end do
end do
end do
! Z2 <-- W6(w,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0) =  &
  W6_(s_a0)%array(i_a0)
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a),&
                     1,&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_x, i_a)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a) * &
                1 * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no2_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no3_x0_type1_eri_c &
  (sx, ix, V2, W7, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sx, ix
real(kind=8), intent(inout) :: V2(*), W7(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sx, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sx

call set_symblock_Xa(sleft, W7, nir, nsym, psym) ! -> Xa (allocate) 
call g_sigma_ccvv_covv_no3_x0_type1_eri_c &
  (sx, ix, h2_i, Xa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xa)

end subroutine g_if_sigma_ccvv_covv_no3_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no3_x0_type1_eri_c &
  (s_x, i_x, V2_, W7_, D2_, nir, nsym, psym, flops)

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
type(symblock1), intent(inout) :: W7_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a3, i_a3, s_a1, i_a1, s_a0, i_a0
! W7(x,a0) += (    1.00000000) V2(x,a2,a3,a1) D2(a3,a1,a2,a0) 
do s_a2 = 0, nir-1
do s_a3 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_a0) == 0 .and. & 
IEOR(s_x,s_a2) == IEOR(s_a3,s_a1) .and. &
IEOR(s_a3,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- D2(a3,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a3, i_a1, i_a2) =  &
  D2_(s_a0, s_a2, s_a1, s_a3)%array(i_a0, i_a2, i_a1, i_a3)
end do
end do
end do
end do
! Z2 <-- V2(x,a2,a3,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a1, i_a2) =  &
  V2_(s_a1, s_a3, s_a2)%array(i_a1, i_a3, i_a2)
end do
end do
end do

! Z3 <-- W7(x,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W7(x,a0)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W7_(s_a0)%array(i_a0) = &
    W7_(s_a0)%array(i_a0) &
  + Z3_(i_a0)
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no3_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no3_x1_type1_eri_c &
  (sb, ib, sx, ix, T2, W7, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sx, ix
real(kind=8), intent(inout) :: T2(*), W7(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sx

call set_symblock_Xa(sleft, W7, nir, nsym, psym) ! -> Xa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no3_x1_type1_eri_c &
  (sb, ib, sx, ix, av2_i, Xa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xa)

end subroutine g_if_sigma_ccvv_covv_no3_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no3_x1_type1_eri_c &
  (s_b, i_b, s_x, i_x, T2_, W7_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W7_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_a, i_a
! S2(w,x,a,b) += (   -2.00000000) T2(w,a0,a,b) W7(x,a0) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_a,s_b) .and. &
IEOR(s_x,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
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
! Z2 <-- W7(x,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0) =  &
  W7_(s_a0)%array(i_a0)
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     1,&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_w, i_a)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                1 * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no3_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no4_x0_type1_eri_c &
  (sx, ix, V2, W8, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sx, ix
real(kind=8), intent(inout) :: V2(*), W8(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sx, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sx

call set_symblock_Xcca(sleft, W8, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_ccvv_covv_no4_x0_type1_eri_c &
  (sx, ix, h2_i, Xcca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccvv_covv_no4_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no4_x0_type1_eri_c &
  (s_x, i_x, V2_, W8_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W8_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_w, i_w, s_a1, i_a1, s_a0, i_a0
! W8(x,w,c0,a0) += (    1.00000000) V2(x,c0,w,a1) D1(a1,a0) 
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_w) == IEOR(s_c0,s_a0) .and. & 
IEOR(s_x,s_c0) == IEOR(s_w,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do
! Z2 <-- V2(x,c0,w,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_c0, i_w) =  &
  V2_(s_a1, s_w, s_c0)%array(i_a1, i_w, i_c0)
end do
end do
end do

! Z3 <-- W8(x,w,c0,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W8(x,w,c0,a0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W8_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w) = &
    W8_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w) &
  + Z3_(i_a0, i_c0, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no4_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no4_x1_type1_eri_c &
  (sb, ib, sx, ix, T2, W8, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sx, ix
real(kind=8), intent(inout) :: T2(*), W8(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sx

call set_symblock_Xcca(sleft, W8, nir, nsym, psym) ! -> Xcca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no4_x1_type1_eri_c &
  (sb, ib, sx, ix, av2_i, Xcca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccvv_covv_no4_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no4_x1_type1_eri_c &
  (s_b, i_b, s_x, i_x, T2_, W8_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_x, s_x
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

integer :: s_a0, i_a0, s_c0, i_c0, s_a, i_a, s_w, i_w
! S2(w,x,a,b) += (    2.00000000) T2(a0,c0,a,b) W8(x,w,c0,a0) 
do s_a0 = 0, nir-1
do s_c0 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_c0) == IEOR(s_a,s_b) .and. &
IEOR(s_x,s_w) == IEOR(s_c0,s_a0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
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
! Z2 <-- W8(x,w,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_c0, i_w) =  &
  W8_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_a, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no4_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no5_x0_type1_eri_c &
  (sx, ix, V2, W9, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sx, ix
real(kind=8), intent(inout) :: V2(*), W9(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sx, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sx

call set_symblock_Xcca(sleft, W9, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_ccvv_covv_no5_x0_type1_eri_c &
  (sx, ix, h2_i, Xcca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccvv_covv_no5_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no5_x0_type1_eri_c &
  (s_x, i_x, V2_, W9_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W9_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_w, i_w, s_c0, i_c0, s_a0, i_a0
! W9(x,w,c0,a0) += (    1.00000000) V2(x,a1,w,c0) D1(a1,a0) 
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_w) == IEOR(s_c0,s_a0) .and. & 
IEOR(s_x,s_a1) == IEOR(s_w,s_c0) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do
! Z2 <-- V2(x,a1,w,c0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w, i_c0) =  &
  V2_(s_c0, s_w, s_a1)%array(i_c0, i_w, i_a1)
end do
end do
end do

! Z3 <-- W9(x,w,c0,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W9(x,w,c0,a0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W9_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w) = &
    W9_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w) &
  + Z3_(i_a0, i_w, i_c0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no5_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no5_x1_type1_eri_c &
  (sb, ib, sx, ix, T2, W9, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sx, ix
real(kind=8), intent(inout) :: T2(*), W9(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sx

call set_symblock_Xcca(sleft, W9, nir, nsym, psym) ! -> Xcca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no5_x1_type1_eri_c &
  (sb, ib, sx, ix, av2_i, Xcca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccvv_covv_no5_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no5_x1_type1_eri_c &
  (s_b, i_b, s_x, i_x, T2_, W9_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_x, s_x
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

integer :: s_a0, i_a0, s_c0, i_c0, s_a, i_a, s_w, i_w
! S2(w,x,a,b) += (   -1.00000000) T2(a0,c0,a,b) W9(x,w,c0,a0) 
do s_a0 = 0, nir-1
do s_c0 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_c0) == IEOR(s_a,s_b) .and. &
IEOR(s_x,s_w) == IEOR(s_c0,s_a0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
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
! Z2 <-- W9(x,w,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_c0, i_w) =  &
  W9_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_a, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no5_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no6_x0_type1_eri_c &
  (sw, iw, V2, W10, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), W10(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_Xa(sleft, W10, nir, nsym, psym) ! -> Xa (allocate) 
call g_sigma_ccvv_covv_no6_x0_type1_eri_c &
  (sw, iw, h2_i, Xa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xa)

end subroutine g_if_sigma_ccvv_covv_no6_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no6_x0_type1_eri_c &
  (s_w, i_w, V2_, W10_, D2_, nir, nsym, psym, flops)

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
type(symblock1), intent(inout) :: W10_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a3, i_a3, s_a1, i_a1, s_a0, i_a0
! W10(w,a0) += (    1.00000000) V2(w,a2,a3,a1) D2(a3,a1,a2,a0) 
do s_a2 = 0, nir-1
do s_a3 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == 0 .and. & 
IEOR(s_w,s_a2) == IEOR(s_a3,s_a1) .and. &
IEOR(s_a3,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- D2(a3,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a3, i_a1, i_a2) =  &
  D2_(s_a0, s_a2, s_a1, s_a3)%array(i_a0, i_a2, i_a1, i_a3)
end do
end do
end do
end do
! Z2 <-- V2(w,a2,a3,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a1, i_a2) =  &
  V2_(s_a1, s_a3, s_a2)%array(i_a1, i_a3, i_a2)
end do
end do
end do

! Z3 <-- W10(w,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W10(w,a0)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W10_(s_a0)%array(i_a0) = &
    W10_(s_a0)%array(i_a0) &
  + Z3_(i_a0)
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no6_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no6_x1_type1_eri_c &
  (sb, ib, sw, iw, T2, W10, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: T2(*), W10(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xa(sleft, W10, nir, nsym, psym) ! -> Xa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no6_x1_type1_eri_c &
  (sb, ib, sw, iw, av2_i, Xa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xa)

end subroutine g_if_sigma_ccvv_covv_no6_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no6_x1_type1_eri_c &
  (s_b, i_b, s_w, i_w, T2_, W10_, S2_, nir, nsym, psym, flops)

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
type(symblock1), intent(inout) :: W10_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_x, i_x, s_a, i_a
! S2(w,x,a,b) += (   -2.00000000) T2(a0,x,a,b) W10(w,a0) 
do s_a0 = 0, nir-1
do s_x = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_x) == IEOR(s_a,s_b) .and. &
IEOR(s_w,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,x,a,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_a, i_a0) =  &
  T2_(s_a, s_x, s_a0)%array(i_a, i_x, i_a0)
end do
end do
end do
! Z2 <-- W10(w,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0) =  &
  W10_(s_a0)%array(i_a0)
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a),&
                     1,&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_x, i_a)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a) * &
                1 * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no6_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no7_x0_type1_eri_c &
  (sx, ix, V2, W11, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sx, ix
real(kind=8), intent(inout) :: V2(*), W11(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sx, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sx

call set_symblock_Xa(sleft, W11, nir, nsym, psym) ! -> Xa (allocate) 
call g_sigma_ccvv_covv_no7_x0_type1_eri_c &
  (sx, ix, h2_i, Xa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xa)

end subroutine g_if_sigma_ccvv_covv_no7_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no7_x0_type1_eri_c &
  (s_x, i_x, V2_, W11_, D2_, nir, nsym, psym, flops)

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
type(symblock1), intent(inout) :: W11_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a3, i_a3, s_a1, i_a1, s_a0, i_a0
! W11(x,a0) += (    1.00000000) V2(x,a2,a3,a1) D2(a3,a1,a2,a0) 
do s_a2 = 0, nir-1
do s_a3 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_a0) == 0 .and. & 
IEOR(s_x,s_a2) == IEOR(s_a3,s_a1) .and. &
IEOR(s_a3,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- D2(a3,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a3, i_a1, i_a2) =  &
  D2_(s_a0, s_a2, s_a1, s_a3)%array(i_a0, i_a2, i_a1, i_a3)
end do
end do
end do
end do
! Z2 <-- V2(x,a2,a3,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a1, i_a2) =  &
  V2_(s_a1, s_a3, s_a2)%array(i_a1, i_a3, i_a2)
end do
end do
end do

! Z3 <-- W11(x,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W11(x,a0)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W11_(s_a0)%array(i_a0) = &
    W11_(s_a0)%array(i_a0) &
  + Z3_(i_a0)
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no7_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no7_x1_type1_eri_c &
  (sb, ib, sx, ix, T2, W11, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sx, ix
real(kind=8), intent(inout) :: T2(*), W11(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sx

call set_symblock_Xa(sleft, W11, nir, nsym, psym) ! -> Xa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no7_x1_type1_eri_c &
  (sb, ib, sx, ix, av2_i, Xa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xa)

end subroutine g_if_sigma_ccvv_covv_no7_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no7_x1_type1_eri_c &
  (s_b, i_b, s_x, i_x, T2_, W11_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W11_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_a, i_a
! S2(w,x,a,b) += (    1.00000000) T2(a0,w,a,b) W11(x,a0) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a,s_b) .and. &
IEOR(s_x,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
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
! Z2 <-- W11(x,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0) =  &
  W11_(s_a0)%array(i_a0)
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     1,&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_w, i_a)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                1 * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no7_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no8_x0_type1_eri_c &
  (sx, ix, V2, W12, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sx, ix
real(kind=8), intent(inout) :: V2(*), W12(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sx, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sx

call set_symblock_Xavv(sleft, W12, nir, nsym, psym) ! -> Xavv (allocate) 
call g_sigma_ccvv_covv_no8_x0_type1_eri_c &
  (sx, ix, h2_i, Xavv, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no8_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no8_x0_type1_eri_c &
  (s_x, i_x, V2_, W12_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W12_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_v0, i_v0, s_a1, i_a1, s_a0, i_a0
! W12(x,a0,v0,a) += (    1.00000000) V2(x,a,v0,a1) D1(a1,a0) 
do s_a = 0, nir-1
do s_v0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_a0) == IEOR(s_v0,s_a) .and. & 
IEOR(s_x,s_a) == IEOR(s_v0,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(x,a,v0,a1) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_v0, i_a1) =  &
  V2_(s_a1, s_v0, s_a)%array(i_a1, i_v0, i_a)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do

! Z3 <-- W12(x,a0,v0,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0))

! W12(x,a0,v0,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W12_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) = &
    W12_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) &
  + Z3_(i_a, i_v0, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no8_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no8_x1_type1_eri_c &
  (sb, ib, sx, ix, T2, W12, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sx, ix
real(kind=8), intent(inout) :: T2(*), W12(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sx

call set_symblock_Xavv(sleft, W12, nir, nsym, psym) ! -> Xavv (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no8_x1_type1_eri_c &
  (sb, ib, sx, ix, av2_i, Xavv, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no8_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no8_x1_type1_eri_c &
  (s_b, i_b, s_x, i_x, T2_, W12_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_x, s_x
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

integer :: s_w, i_w, s_a0, i_a0, s_v0, i_v0, s_a, i_a
! S2(w,x,a,b) += (    1.00000000) T2(w,a0,v0,b) W12(x,a0,v0,a) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_v0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_v0,s_b) .and. &
IEOR(s_x,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W12(x,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0, i_v0) =  &
  W12_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0)
end do
end do
end do
! Z2 <-- T2(w,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_v0, i_w) =  &
  T2_(s_v0, s_a0, s_w)%array(i_v0, i_a0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_a, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no8_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no9_x0_type1_eri_c &
  (sw, iw, V2, W13, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), W13(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_Xavv(sleft, W13, nir, nsym, psym) ! -> Xavv (allocate) 
call g_sigma_ccvv_covv_no9_x0_type1_eri_c &
  (sw, iw, h2_i, Xavv, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no9_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no9_x0_type1_eri_c &
  (s_w, i_w, V2_, W13_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W13_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_v0, i_v0, s_a1, i_a1, s_a0, i_a0
! W13(w,a0,v0,a) += (    1.00000000) V2(w,a,v0,a1) D1(a1,a0) 
do s_a = 0, nir-1
do s_v0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_v0,s_a) .and. & 
IEOR(s_w,s_a) == IEOR(s_v0,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(w,a,v0,a1) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_v0, i_a1) =  &
  V2_(s_a1, s_v0, s_a)%array(i_a1, i_v0, i_a)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do

! Z3 <-- W13(w,a0,v0,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0))

! W13(w,a0,v0,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W13_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) = &
    W13_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) &
  + Z3_(i_a, i_v0, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no9_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no9_x1_type1_eri_c &
  (sb, ib, sw, iw, T2, W13, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: T2(*), W13(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xavv(sleft, W13, nir, nsym, psym) ! -> Xavv (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no9_x1_type1_eri_c &
  (sb, ib, sw, iw, av2_i, Xavv, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no9_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no9_x1_type1_eri_c &
  (s_b, i_b, s_w, i_w, T2_, W13_, S2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W13_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_a0, i_a0, s_v0, i_v0, s_a, i_a
! S2(w,x,a,b) += (   -2.00000000) T2(x,a0,v0,b) W13(w,a0,v0,a) 
do s_x = 0, nir-1
do s_a0 = 0, nir-1
do s_v0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_x,s_a0) == IEOR(s_v0,s_b) .and. &
IEOR(s_w,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W13(w,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0, i_v0) =  &
  W13_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0)
end do
end do
end do
! Z2 <-- T2(x,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_v0, i_x) =  &
  T2_(s_v0, s_a0, s_x)%array(i_v0, i_a0, i_x)
end do
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_a, i_x)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no9_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no10_x0_type1_eri_c &
  (sw, iw, V2, W18, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), W18(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_Xavv(sleft, W18, nir, nsym, psym) ! -> Xavv (allocate) 
call g_sigma_ccvv_covv_no10_x0_type1_eri_c &
  (sw, iw, h2_i, Xavv, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no10_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no10_x0_type1_eri_c &
  (s_w, i_w, V2_, W18_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W18_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_v0, i_v0, s_a, i_a, s_a0, i_a0
! W18(w,a0,v0,a) += (    1.00000000) V2(w,a1,v0,a) D1(a1,a0) 
do s_a1 = 0, nir-1
do s_v0 = 0, nir-1
do s_a = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_v0,s_a) .and. & 
IEOR(s_w,s_a1) == IEOR(s_v0,s_a) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(w,a1,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z1_(i_v0, i_a, i_a1) =  &
  V2_(s_a, s_v0, s_a1)%array(i_a, i_v0, i_a1)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do

! Z3 <-- W18(w,a0,v0,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a))

! W18(w,a0,v0,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W18_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) = &
    W18_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) &
  + Z3_(i_v0, i_a, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no10_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no10_x1_type1_eri_c &
  (sb, ib, sw, iw, T2, W18, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: T2(*), W18(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xavv(sleft, W18, nir, nsym, psym) ! -> Xavv (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no10_x1_type1_eri_c &
  (sb, ib, sw, iw, av2_i, Xavv, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no10_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no10_x1_type1_eri_c &
  (s_b, i_b, s_w, i_w, T2_, W18_, S2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W18_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_a0, i_a0, s_v0, i_v0, s_a, i_a
! S2(w,x,a,b) += (    1.00000000) T2(x,a0,v0,b) W18(w,a0,v0,a) 
do s_x = 0, nir-1
do s_a0 = 0, nir-1
do s_v0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_x,s_a0) == IEOR(s_v0,s_b) .and. &
IEOR(s_w,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W18(w,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0, i_v0) =  &
  W18_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0)
end do
end do
end do
! Z2 <-- T2(x,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_v0, i_x) =  &
  T2_(s_v0, s_a0, s_x)%array(i_v0, i_a0, i_x)
end do
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_a, i_x)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no10_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no11_x0_type1_eri_c &
  (sx, ix, V2, W19, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sx, ix
real(kind=8), intent(inout) :: V2(*), W19(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sx, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sx

call set_symblock_Xavv(sleft, W19, nir, nsym, psym) ! -> Xavv (allocate) 
call g_sigma_ccvv_covv_no11_x0_type1_eri_c &
  (sx, ix, h2_i, Xavv, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no11_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no11_x0_type1_eri_c &
  (s_x, i_x, V2_, W19_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W19_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_v0, i_v0, s_a, i_a, s_a0, i_a0
! W19(x,a0,v0,a) += (    1.00000000) V2(x,a1,v0,a) D1(a1,a0) 
do s_a1 = 0, nir-1
do s_v0 = 0, nir-1
do s_a = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_a0) == IEOR(s_v0,s_a) .and. & 
IEOR(s_x,s_a1) == IEOR(s_v0,s_a) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(x,a1,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z1_(i_v0, i_a, i_a1) =  &
  V2_(s_a, s_v0, s_a1)%array(i_a, i_v0, i_a1)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do

! Z3 <-- W19(x,a0,v0,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a))

! W19(x,a0,v0,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W19_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) = &
    W19_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) &
  + Z3_(i_v0, i_a, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no11_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no11_x1_type1_eri_c &
  (sb, ib, sx, ix, T2, W19, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sx, ix
real(kind=8), intent(inout) :: T2(*), W19(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sx

call set_symblock_Xavv(sleft, W19, nir, nsym, psym) ! -> Xavv (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no11_x1_type1_eri_c &
  (sb, ib, sx, ix, av2_i, Xavv, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no11_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no11_x1_type1_eri_c &
  (s_b, i_b, s_x, i_x, T2_, W19_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W19_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_v0, i_v0, s_a, i_a
! S2(w,x,a,b) += (   -2.00000000) T2(w,a0,v0,b) W19(x,a0,v0,a) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_v0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_v0,s_b) .and. &
IEOR(s_x,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W19(x,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0, i_v0) =  &
  W19_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0)
end do
end do
end do
! Z2 <-- T2(w,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_v0, i_w) =  &
  T2_(s_v0, s_a0, s_w)%array(i_v0, i_a0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_a, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no11_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no12_x0_type1_eri_c &
  (sw, iw, V2, W20, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), W20(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_Xavv(sleft, W20, nir, nsym, psym) ! -> Xavv (allocate) 
call g_sigma_ccvv_covv_no12_x0_type1_eri_c &
  (sw, iw, h2_i, Xavv, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no12_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no12_x0_type1_eri_c &
  (s_w, i_w, V2_, W20_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W20_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_v0, i_v0, s_a1, i_a1, s_a0, i_a0
! W20(w,a0,v0,a) += (    1.00000000) V2(w,a,v0,a1) D1(a1,a0) 
do s_a = 0, nir-1
do s_v0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_v0,s_a) .and. & 
IEOR(s_w,s_a) == IEOR(s_v0,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(w,a,v0,a1) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_v0, i_a1) =  &
  V2_(s_a1, s_v0, s_a)%array(i_a1, i_v0, i_a)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do

! Z3 <-- W20(w,a0,v0,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0))

! W20(w,a0,v0,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W20_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) = &
    W20_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) &
  + Z3_(i_a, i_v0, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no12_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no12_x1_type1_eri_c &
  (sb, ib, sw, iw, T2, W20, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: T2(*), W20(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xavv(sleft, W20, nir, nsym, psym) ! -> Xavv (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no12_x1_type1_eri_c &
  (sb, ib, sw, iw, av2_i, Xavv, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no12_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no12_x1_type1_eri_c &
  (s_b, i_b, s_w, i_w, T2_, W20_, S2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W20_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_x, i_x, s_v0, i_v0, s_a, i_a
! S2(w,x,a,b) += (    4.00000000) T2(a0,x,v0,b) W20(w,a0,v0,a) 
do s_a0 = 0, nir-1
do s_x = 0, nir-1
do s_v0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_x) == IEOR(s_v0,s_b) .and. &
IEOR(s_w,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W20(w,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0, i_v0) =  &
  W20_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0)
end do
end do
end do
! Z2 <-- T2(a0,x,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_v0, i_x) =  &
  T2_(s_v0, s_x, s_a0)%array(i_v0, i_x, i_a0)
end do
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_a, i_x)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no12_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no13_x0_type1_eri_c &
  (sx, ix, V2, W21, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sx, ix
real(kind=8), intent(inout) :: V2(*), W21(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sx, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sx

call set_symblock_Xavv(sleft, W21, nir, nsym, psym) ! -> Xavv (allocate) 
call g_sigma_ccvv_covv_no13_x0_type1_eri_c &
  (sx, ix, h2_i, Xavv, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no13_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no13_x0_type1_eri_c &
  (s_x, i_x, V2_, W21_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W21_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_v0, i_v0, s_a1, i_a1, s_a0, i_a0
! W21(x,a0,v0,a) += (    1.00000000) V2(x,a,v0,a1) D1(a1,a0) 
do s_a = 0, nir-1
do s_v0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_a0) == IEOR(s_v0,s_a) .and. & 
IEOR(s_x,s_a) == IEOR(s_v0,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(x,a,v0,a1) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_v0, i_a1) =  &
  V2_(s_a1, s_v0, s_a)%array(i_a1, i_v0, i_a)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do

! Z3 <-- W21(x,a0,v0,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0))

! W21(x,a0,v0,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W21_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) = &
    W21_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) &
  + Z3_(i_a, i_v0, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a)*psym(I_LENGTH,I_V, s_v0) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no13_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no13_x1_type1_eri_c &
  (sb, ib, sx, ix, T2, W21, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sx, ix
real(kind=8), intent(inout) :: T2(*), W21(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sx

call set_symblock_Xavv(sleft, W21, nir, nsym, psym) ! -> Xavv (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no13_x1_type1_eri_c &
  (sb, ib, sx, ix, av2_i, Xavv, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no13_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no13_x1_type1_eri_c &
  (s_b, i_b, s_x, i_x, T2_, W21_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_x, s_x
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
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_v0, i_v0, s_a, i_a
! S2(w,x,a,b) += (   -2.00000000) T2(a0,w,v0,b) W21(x,a0,v0,a) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_v0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_w) == IEOR(s_v0,s_b) .and. &
IEOR(s_x,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W21(x,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0, i_v0) =  &
  W21_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0)
end do
end do
end do
! Z2 <-- T2(a0,w,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_v0, i_w) =  &
  T2_(s_v0, s_w, s_a0)%array(i_v0, i_w, i_a0)
end do
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_a, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no13_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no14_x0_type1_eri_c &
  (sw, iw, V2, W26, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), W26(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_Xavv(sleft, W26, nir, nsym, psym) ! -> Xavv (allocate) 
call g_sigma_ccvv_covv_no14_x0_type1_eri_c &
  (sw, iw, h2_i, Xavv, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no14_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no14_x0_type1_eri_c &
  (s_w, i_w, V2_, W26_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W26_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_v0, i_v0, s_a, i_a, s_a0, i_a0
! W26(w,a0,v0,a) += (    1.00000000) V2(w,a1,v0,a) D1(a1,a0) 
do s_a1 = 0, nir-1
do s_v0 = 0, nir-1
do s_a = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_v0,s_a) .and. & 
IEOR(s_w,s_a1) == IEOR(s_v0,s_a) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(w,a1,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z1_(i_v0, i_a, i_a1) =  &
  V2_(s_a, s_v0, s_a1)%array(i_a, i_v0, i_a1)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do

! Z3 <-- W26(w,a0,v0,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a))

! W26(w,a0,v0,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W26_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) = &
    W26_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) &
  + Z3_(i_v0, i_a, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no14_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no14_x1_type1_eri_c &
  (sb, ib, sw, iw, T2, W26, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sw, iw
real(kind=8), intent(inout) :: T2(*), W26(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xavv(sleft, W26, nir, nsym, psym) ! -> Xavv (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no14_x1_type1_eri_c &
  (sb, ib, sw, iw, av2_i, Xavv, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no14_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no14_x1_type1_eri_c &
  (s_b, i_b, s_w, i_w, T2_, W26_, S2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W26_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_x, i_x, s_v0, i_v0, s_a, i_a
! S2(w,x,a,b) += (   -2.00000000) T2(a0,x,v0,b) W26(w,a0,v0,a) 
do s_a0 = 0, nir-1
do s_x = 0, nir-1
do s_v0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_x) == IEOR(s_v0,s_b) .and. &
IEOR(s_w,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W26(w,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0, i_v0) =  &
  W26_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0)
end do
end do
end do
! Z2 <-- T2(a0,x,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_v0, i_x) =  &
  T2_(s_v0, s_x, s_a0)%array(i_v0, i_x, i_a0)
end do
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_a, i_x)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no14_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no15_x0_type1_eri_c &
  (sx, ix, V2, W27, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sx, ix
real(kind=8), intent(inout) :: V2(*), W27(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sx, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sx

call set_symblock_Xavv(sleft, W27, nir, nsym, psym) ! -> Xavv (allocate) 
call g_sigma_ccvv_covv_no15_x0_type1_eri_c &
  (sx, ix, h2_i, Xavv, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no15_x0_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no15_x0_type1_eri_c &
  (s_x, i_x, V2_, W27_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W27_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_v0, i_v0, s_a, i_a, s_a0, i_a0
! W27(x,a0,v0,a) += (    1.00000000) V2(x,a1,v0,a) D1(a1,a0) 
do s_a1 = 0, nir-1
do s_v0 = 0, nir-1
do s_a = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_a0) == IEOR(s_v0,s_a) .and. & 
IEOR(s_x,s_a1) == IEOR(s_v0,s_a) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(x,a1,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z1_(i_v0, i_a, i_a1) =  &
  V2_(s_a, s_v0, s_a1)%array(i_a, i_v0, i_a1)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do

! Z3 <-- W27(x,a0,v0,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a))

! W27(x,a0,v0,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W27_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) = &
    W27_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0) &
  + Z3_(i_v0, i_a, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no15_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no15_x1_type1_eri_c &
  (sb, ib, sx, ix, T2, W27, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sx, ix
real(kind=8), intent(inout) :: T2(*), W27(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sx

call set_symblock_Xavv(sleft, W27, nir, nsym, psym) ! -> Xavv (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no15_x1_type1_eri_c &
  (sb, ib, sx, ix, av2_i, Xavv, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xavv)

end subroutine g_if_sigma_ccvv_covv_no15_x1_type1_eri_c



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
subroutine g_sigma_ccvv_covv_no15_x1_type1_eri_c &
  (s_b, i_b, s_x, i_x, T2_, W27_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_x, s_x
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W27_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_v0, i_v0, s_a, i_a
! S2(w,x,a,b) += (    1.00000000) T2(a0,w,v0,b) W27(x,a0,v0,a) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_v0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_w) == IEOR(s_v0,s_b) .and. &
IEOR(s_x,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W27(x,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0, i_v0) =  &
  W27_(s_a, s_v0, s_a0)%array(i_a, i_v0, i_a0)
end do
end do
end do
! Z2 <-- T2(a0,w,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_v0, i_w) =  &
  T2_(s_v0, s_w, s_a0)%array(i_v0, i_w, i_a0)
end do
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_a, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no15_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no0_x0_type1_eri_v &
  (sb, ib, V2, W14, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: V2(*), W14(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sb

call set_symblock_Xcav(sleft, W14, nir, nsym, psym) ! -> Xcav (allocate) 
call g_sigma_ccvv_covv_no0_x0_type1_eri_v &
  (sb, ib, h2_i, Xcav, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcav)

end subroutine g_if_sigma_ccvv_covv_no0_x0_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no0_x0_type1_eri_v &
  (s_b, i_b, V2_, W14_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W14_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_v0, i_v0, s_a1, i_a1, s_a0, i_a0
! W14(x,a0,v0,b) += (    1.00000000) V2(b,x,v0,a1) D1(a1,a0) 
do s_x = 0, nir-1
do s_v0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_a0) == IEOR(s_v0,s_b) .and. & 
IEOR(s_b,s_x) == IEOR(s_v0,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(b,x,v0,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_v0, i_a1) =  &
  V2_(s_a1, s_v0, s_x)%array(i_a1, i_v0, i_x)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do

! Z3 <-- W14(x,a0,v0,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_v0),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_v0))

! W14(x,a0,v0,b)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
W14_(s_v0, s_a0, s_x)%array(i_v0, i_a0, i_x) = &
    W14_(s_v0, s_a0, s_x)%array(i_v0, i_a0, i_x) &
  + Z3_(i_x, i_v0, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_v0) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no0_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no0_x1_type1_eri_v &
  (sa, ia, sb, ib, T2, W14, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sb, ib
real(kind=8), intent(inout) :: T2(*), W14(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xcav(sleft, W14, nir, nsym, psym) ! -> Xcav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no0_x1_type1_eri_v &
  (sa, ia, sb, ib, av2_i, Xcav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcav)

end subroutine g_if_sigma_ccvv_covv_no0_x1_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no0_x1_type1_eri_v &
  (s_a, i_a, s_b, i_b, T2_, W14_, S2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W14_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_v0, i_v0, s_x, i_x
! S2(w,x,a,b) += (   -2.00000000) T2(w,a0,v0,a) W14(x,a0,v0,b) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_v0 = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_v0,s_a) .and. &
IEOR(s_x,s_a0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W14(x,a0,v0,b) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_a0, i_v0) =  &
  W14_(s_v0, s_a0, s_x)%array(i_v0, i_a0, i_x)
end do
end do
end do
! Z2 <-- T2(w,a0,v0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_v0, i_w) =  &
  T2_(s_v0, s_a0, s_w)%array(i_v0, i_a0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_x, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no0_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no1_x0_type1_eri_v &
  (sb, ib, V2, W15, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: V2(*), W15(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sb

call set_symblock_Xcav(sleft, W15, nir, nsym, psym) ! -> Xcav (allocate) 
call g_sigma_ccvv_covv_no1_x0_type1_eri_v &
  (sb, ib, h2_i, Xcav, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcav)

end subroutine g_if_sigma_ccvv_covv_no1_x0_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no1_x0_type1_eri_v &
  (s_b, i_b, V2_, W15_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W15_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_v0, i_v0, s_a1, i_a1, s_a0, i_a0
! W15(w,a0,v0,b) += (    1.00000000) V2(b,w,v0,a1) D1(a1,a0) 
do s_w = 0, nir-1
do s_v0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_v0,s_b) .and. & 
IEOR(s_b,s_w) == IEOR(s_v0,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(b,w,v0,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_v0, i_a1) =  &
  V2_(s_a1, s_v0, s_w)%array(i_a1, i_v0, i_w)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do

! Z3 <-- W15(w,a0,v0,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0))

! W15(w,a0,v0,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
W15_(s_v0, s_a0, s_w)%array(i_v0, i_a0, i_w) = &
    W15_(s_v0, s_a0, s_w)%array(i_v0, i_a0, i_w) &
  + Z3_(i_w, i_v0, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no1_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no1_x1_type1_eri_v &
  (sa, ia, sb, ib, T2, W15, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sb, ib
real(kind=8), intent(inout) :: T2(*), W15(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xcav(sleft, W15, nir, nsym, psym) ! -> Xcav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no1_x1_type1_eri_v &
  (sa, ia, sb, ib, av2_i, Xcav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcav)

end subroutine g_if_sigma_ccvv_covv_no1_x1_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no1_x1_type1_eri_v &
  (s_a, i_a, s_b, i_b, T2_, W15_, S2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W15_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_a0, i_a0, s_v0, i_v0, s_w, i_w
! S2(w,x,a,b) += (    1.00000000) T2(x,a0,v0,a) W15(w,a0,v0,b) 
do s_x = 0, nir-1
do s_a0 = 0, nir-1
do s_v0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_x,s_a0) == IEOR(s_v0,s_a) .and. &
IEOR(s_w,s_a0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- T2(x,a0,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_a0, i_v0) =  &
  T2_(s_v0, s_a0, s_x)%array(i_v0, i_a0, i_x)
end do
end do
end do
! Z2 <-- W15(w,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_v0, i_w) =  &
  W15_(s_v0, s_a0, s_w)%array(i_v0, i_a0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_x, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no1_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no2_x0_type1_eri_v &
  (sb, ib, sv0, iv0, V2, W16, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W16(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xca(sleft, W16, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_ccvv_covv_no2_x0_type1_eri_v &
  (sb, ib, sv0, iv0, h2_i, Xca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no2_x0_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no2_x0_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, V2_, W16_, D1_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W16_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a1, i_a1, s_a0, i_a0
! W16(w,a0,v0,b) += (    1.00000000) V2(v0,b,w,a1) D1(a1,a0) 
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_v0,s_b) .and. & 
IEOR(s_v0,s_b) == IEOR(s_w,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do
! Z2 <-- V2(v0,b,w,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w) =  &
  V2_(s_a1, s_w, s_b)%array(i_a1, i_w, i_b)
end do
end do

! Z3 <-- W16(w,a0,v0,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W16(w,a0,v0,b)  <-- Z3
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
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no2_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no2_x1_type1_eri_v &
  (sb, ib, sv0, iv0, T2, W16, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W16(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xca(sleft, W16, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no2_x1_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no2_x1_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no2_x1_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, W16_, S2_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W16_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_x, i_x, s_a, i_a, s_w, i_w
! S2(w,x,a,b) += (   -2.00000000) T2(a0,x,a,v0) W16(w,a0,v0,b) 
do s_a0 = 0, nir-1
do s_x = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_x) == IEOR(s_a,s_v0) .and. &
IEOR(s_w,s_a0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,x,a,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_a, i_a0) =  &
  T2_(s_a, s_x, s_a0)%array(i_a, i_x, i_a0)
end do
end do
end do
! Z2 <-- W16(w,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  W16_(s_a0, s_w)%array(i_a0, i_w)
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_x, i_a, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a) * &
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

end subroutine g_sigma_ccvv_covv_no2_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no3_x0_type1_eri_v &
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

call set_symblock_Xca(sleft, W17, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_ccvv_covv_no3_x0_type1_eri_v &
  (sb, ib, sv0, iv0, h2_i, Xca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no3_x0_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no3_x0_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, V2_, W17_, D1_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W17_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_a1, i_a1, s_a0, i_a0
! W17(x,a0,v0,b) += (    1.00000000) V2(v0,b,x,a1) D1(a1,a0) 
do s_x = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_a0) == IEOR(s_v0,s_b) .and. & 
IEOR(s_v0,s_b) == IEOR(s_x,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do
! Z2 <-- V2(v0,b,x,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_x) =  &
  V2_(s_a1, s_x, s_b)%array(i_a1, i_x, i_b)
end do
end do

! Z3 <-- W17(x,a0,v0,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W17(x,a0,v0,b)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W17_(s_a0, s_x)%array(i_a0, i_x) = &
    W17_(s_a0, s_x)%array(i_a0, i_x) &
  + Z3_(i_a0, i_x)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no3_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no3_x1_type1_eri_v &
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

call set_symblock_Xca(sleft, W17, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no3_x1_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no3_x1_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no3_x1_type1_eri_v &
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

integer :: s_a0, i_a0, s_w, i_w, s_a, i_a, s_x, i_x
! S2(w,x,a,b) += (    1.00000000) T2(a0,w,a,v0) W17(x,a0,v0,b) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a,s_v0) .and. &
IEOR(s_x,s_a0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
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
! Z2 <-- W17(x,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x) =  &
  W17_(s_a0, s_x)%array(i_a0, i_x)
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_w, i_a, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no3_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no4_x0_type1_eri_v &
  (sb, ib, sv0, iv0, V2, W22, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W22(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xca(sleft, W22, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_ccvv_covv_no4_x0_type1_eri_v &
  (sb, ib, sv0, iv0, h2_i, Xca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no4_x0_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no4_x0_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, V2_, W22_, D1_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W22_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_a1, i_a1, s_a0, i_a0
! W22(x,a0,v0,b) += (    1.00000000) V2(b,x,v0,a1) D1(a1,a0) 
do s_x = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_a0) == IEOR(s_v0,s_b) .and. & 
IEOR(s_b,s_x) == IEOR(s_v0,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do
! Z2 <-- V2(b,x,v0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_x) =  &
  V2_(s_a1, s_v0, s_x)%array(i_a1, i_v0, i_x)
end do
end do

! Z3 <-- W22(x,a0,v0,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W22(x,a0,v0,b)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W22_(s_a0, s_x)%array(i_a0, i_x) = &
    W22_(s_a0, s_x)%array(i_a0, i_x) &
  + Z3_(i_a0, i_x)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no4_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no4_x1_type1_eri_v &
  (sb, ib, sv0, iv0, T2, W22, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W22(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xca(sleft, W22, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no4_x1_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no4_x1_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no4_x1_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, W22_, S2_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W22_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_a, i_a, s_x, i_x
! S2(w,x,a,b) += (    4.00000000) T2(w,a0,a,v0) W22(x,a0,v0,b) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_a,s_v0) .and. &
IEOR(s_x,s_a0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
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
! Z2 <-- W22(x,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x) =  &
  W22_(s_a0, s_x)%array(i_a0, i_x)
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_w, i_a, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no4_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no5_x0_type1_eri_v &
  (sb, ib, sv0, iv0, V2, W23, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W23(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xca(sleft, W23, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_ccvv_covv_no5_x0_type1_eri_v &
  (sb, ib, sv0, iv0, h2_i, Xca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no5_x0_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no5_x0_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, V2_, W23_, D1_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W23_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a1, i_a1, s_a0, i_a0
! W23(w,a0,v0,b) += (    1.00000000) V2(b,w,v0,a1) D1(a1,a0) 
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_v0,s_b) .and. & 
IEOR(s_b,s_w) == IEOR(s_v0,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do
! Z2 <-- V2(b,w,v0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w) =  &
  V2_(s_a1, s_v0, s_w)%array(i_a1, i_v0, i_w)
end do
end do

! Z3 <-- W23(w,a0,v0,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W23(w,a0,v0,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W23_(s_a0, s_w)%array(i_a0, i_w) = &
    W23_(s_a0, s_w)%array(i_a0, i_w) &
  + Z3_(i_a0, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no5_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no5_x1_type1_eri_v &
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
sleft = IEOR(sv0,sb)

call set_symblock_Xca(sleft, W23, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no5_x1_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no5_x1_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no5_x1_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, W23_, S2_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W23_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_a0, i_a0, s_a, i_a, s_w, i_w
! S2(w,x,a,b) += (   -2.00000000) T2(x,a0,a,v0) W23(w,a0,v0,b) 
do s_x = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_x,s_a0) == IEOR(s_a,s_v0) .and. &
IEOR(s_w,s_a0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(x,a0,a,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_a, i_a0) =  &
  T2_(s_a, s_a0, s_x)%array(i_a, i_a0, i_x)
end do
end do
end do
! Z2 <-- W23(w,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  W23_(s_a0, s_w)%array(i_a0, i_w)
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_x, i_a, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a) * &
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

end subroutine g_sigma_ccvv_covv_no5_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no6_x0_type1_eri_v &
  (sb, ib, sv0, iv0, V2, W24, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W24(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xca(sleft, W24, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_ccvv_covv_no6_x0_type1_eri_v &
  (sb, ib, sv0, iv0, h2_i, Xca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no6_x0_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no6_x0_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, V2_, W24_, D1_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W24_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a1, i_a1, s_a0, i_a0
! W24(w,a0,v0,b) += (    1.00000000) V2(v0,b,w,a1) D1(a1,a0) 
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_v0,s_b) .and. & 
IEOR(s_v0,s_b) == IEOR(s_w,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do
! Z2 <-- V2(v0,b,w,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w) =  &
  V2_(s_a1, s_w, s_b)%array(i_a1, i_w, i_b)
end do
end do

! Z3 <-- W24(w,a0,v0,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W24(w,a0,v0,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W24_(s_a0, s_w)%array(i_a0, i_w) = &
    W24_(s_a0, s_w)%array(i_a0, i_w) &
  + Z3_(i_a0, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no6_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no6_x1_type1_eri_v &
  (sb, ib, sv0, iv0, T2, W24, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W24(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sv0,sb)

call set_symblock_Xca(sleft, W24, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no6_x1_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no6_x1_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no6_x1_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, W24_, S2_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W24_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_a0, i_a0, s_a, i_a, s_w, i_w
! S2(w,x,a,b) += (    1.00000000) T2(x,a0,a,v0) W24(w,a0,v0,b) 
do s_x = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_x,s_a0) == IEOR(s_a,s_v0) .and. &
IEOR(s_w,s_a0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(x,a0,a,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_a, i_a0) =  &
  T2_(s_a, s_a0, s_x)%array(i_a, i_a0, i_x)
end do
end do
end do
! Z2 <-- W24(w,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  W24_(s_a0, s_w)%array(i_a0, i_w)
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_x, i_a, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_V, s_a) * &
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

end subroutine g_sigma_ccvv_covv_no6_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no7_x0_type1_eri_v &
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

call set_symblock_Xca(sleft, W25, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_ccvv_covv_no7_x0_type1_eri_v &
  (sb, ib, sv0, iv0, h2_i, Xca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no7_x0_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no7_x0_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, V2_, W25_, D1_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W25_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_a1, i_a1, s_a0, i_a0
! W25(x,a0,v0,b) += (    1.00000000) V2(v0,b,x,a1) D1(a1,a0) 
do s_x = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_a0) == IEOR(s_v0,s_b) .and. & 
IEOR(s_v0,s_b) == IEOR(s_x,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
end do
! Z2 <-- V2(v0,b,x,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_x) =  &
  V2_(s_a1, s_x, s_b)%array(i_a1, i_x, i_b)
end do
end do

! Z3 <-- W25(x,a0,v0,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W25(x,a0,v0,b)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W25_(s_a0, s_x)%array(i_a0, i_x) = &
    W25_(s_a0, s_x)%array(i_a0, i_x) &
  + Z3_(i_a0, i_x)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no7_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccvv_covv_no7_x1_type1_eri_v &
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

call set_symblock_Xca(sleft, W25, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccvv_covv_no7_x1_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_ccvv_covv_no7_x1_type1_eri_v



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
subroutine g_sigma_ccvv_covv_no7_x1_type1_eri_v &
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

integer :: s_w, i_w, s_a0, i_a0, s_a, i_a, s_x, i_x
! S2(w,x,a,b) += (   -2.00000000) T2(w,a0,a,v0) W25(x,a0,v0,b) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a,s_b) .and. & 
IEOR(s_w,s_a0) == IEOR(s_a,s_v0) .and. &
IEOR(s_x,s_a0) == IEOR(s_v0,s_b)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
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
! Z2 <-- W25(x,a0,v0,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x) =  &
  W25_(s_a0, s_x)%array(i_a0, i_x)
end do
end do

! Z3 <-- S2(w,x,a,b) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a))

! S2(w,x,a,b)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) = &
    S2_(s_a, s_x, s_w)%array(i_a, i_x, i_w) &
  + Z3_(i_w, i_a, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccvv_covv_no7_x1_type1_eri_v

