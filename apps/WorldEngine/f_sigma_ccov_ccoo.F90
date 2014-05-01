#include <sci/icmr/fsrc/f_mr.fh>


!      _/_/_/_/                            _/             
!     _/        _/_/    _/_/_/  _/_/    _/_/_/_/    _/_/  
!    _/_/_/  _/_/_/_/  _/    _/    _/    _/      _/    _/ 
!   _/      _/        _/    _/    _/    _/      _/    _/  
!  _/        _/_/_/  _/    _/    _/      _/_/    _/_/     

!                                    Generated date : Sun Apr 20 10:26:24 2014



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no0_x0_type0_noeri &
  (sa1, ia1, T2, W0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa1, ia1
real(kind=8), intent(inout) :: T2(*), W0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xccaa(sleft, W0, nir, nsym, psym) ! -> Xccaa (allocate) 
call g_sigma_ccov_ccoo_no0_x0_type0_noeri &
  (sa1, ia1, av2_i, Xccaa, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xccaa)

end subroutine g_if_sigma_ccov_ccoo_no0_x0_type0_noeri



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
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ccoo_no0_x0_type0_noeri &
  (s_a1, i_a1, T2_, W0_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a0, i_a0, s_i, i_i, s_a2, i_a2
! W0(x,w,i,a2) += (    1.00000000) T2(x,w,a0,a1) D2(i,a1,a2,a0) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_a2 = 0, nir-1
if( &
IEOR(s_x,s_w) == IEOR(s_i,s_a2) .and. & 
IEOR(s_x,s_w) == IEOR(s_a0,s_a1) .and. &
IEOR(s_i,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(i,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a2, i_a0) =  &
  D2_(s_a0, s_a2, s_a1, s_i)%array(i_a0, i_a2, i_a1, i_i)
end do
end do
end do
! Z2 <-- T2(x,w,a0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x, i_w) =  &
  T2_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x)
end do
end do
end do

! Z3 <-- W0(x,w,i,a2) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2))

! W0(x,w,i,a2)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W0_(s_a2, s_i, s_w, s_x)%array(i_a2, i_i, i_w, i_x) = &
    W0_(s_a2, s_i, s_w, s_x)%array(i_a2, i_i, i_w, i_x) &
  + Z3_(i_i, i_a2, i_x, i_w)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) * &
                psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no0_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no0_x1_type0_noeri &
  (sa, ia, W0, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: W0(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xccaa(sleft, W0, nir, nsym, psym) ! -> Xccaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no0_x1_type0_noeri &
  (sa, ia, Xccaa, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xccaa)

end subroutine g_if_sigma_ccov_ccoo_no0_x1_type0_noeri



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
subroutine g_sigma_ccov_ccoo_no0_x1_type0_noeri &
  (s_a, i_a, W0_, S2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_x, i_x, s_w, i_w, s_i, i_i
! S2(w,x,i,a) += (    2.00000000) Fc1(a2,a) W0(x,w,i,a2) 
do s_a2 = 0, nir-1
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a2,s_a) == 0 .and. &
IEOR(s_x,s_w) == IEOR(s_i,s_a2)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- W0(x,w,i,a2) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_w, i_i, i_a2) =  &
  W0_(s_a2, s_i, s_w, s_x)%array(i_a2, i_i, i_w, i_x)
end do
end do
end do
end do
! Z2 <-- Fc1(a2,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2) =  &
  Fc1_(s_a, s_a2)%array(i_a, i_a2)
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_O, s_a2),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
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
                psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no0_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no1_x0_type0_noeri &
  (sa, ia, W1, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: W1(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa

call set_symblock_Xa(sleft, W1, nir, nsym, psym) ! -> Xa (allocate) 
call g_sigma_ccov_ccoo_no1_x0_type0_noeri &
  (sa, ia, Xa, d1, fc1, nir, nsym, psym, flops)

deallocate(Xa)

end subroutine g_if_sigma_ccov_ccoo_no1_x0_type0_noeri



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
subroutine g_sigma_ccov_ccoo_no1_x0_type0_noeri &
  (s_a, i_a, W1_, D1_, Fc1_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W1_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a0, i_a0
! W1(a0,a) += (    1.00000000) D1(a1,a0) Fc1(a1,a) 
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_a0,s_a) == 0 .and. & 
IEOR(s_a1,s_a0) == 0 .and. &
IEOR(s_a1,s_a) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
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
! Z2 <-- Fc1(a1,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1) =  &
  Fc1_(s_a, s_a1)%array(i_a, i_a1)
end do

! Z3 <-- W1(a0,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W1(a0,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W1_(s_a0)%array(i_a0) = &
    W1_(s_a0)%array(i_a0) &
  + Z3_(i_a0)
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no1_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no1_x1_type0_noeri &
  (sa, ia, si, ii, T2, W1, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, si, ii
real(kind=8), intent(inout) :: T2(*), W1(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(si, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xa(sleft, W1, nir, nsym, psym) ! -> Xa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no1_x1_type0_noeri &
  (sa, ia, si, ii, av2_i, Xa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xa)

end subroutine g_if_sigma_ccov_ccoo_no1_x1_type0_noeri



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
subroutine g_sigma_ccov_ccoo_no1_x1_type0_noeri &
  (s_a, i_a, s_i, i_i, T2_, W1_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W1_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a0, i_a0
! S2(w,x,i,a) += (   -4.00000000) T2(x,w,a0,i) W1(a0,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_a0,s_i) .and. &
IEOR(s_a0,s_a) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(x,w,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_w, i_a0) =  &
  T2_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x)
end do
end do
end do
! Z2 <-- W1(a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0) =  &
  W1_(s_a0)%array(i_a0)
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w),&
                     1,&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) * &
                1 * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no1_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no2_x0_type0_noeri &
  (sa, ia, sa0, ia0, W2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: W2
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa0,sa)

call g_sigma_ccov_ccoo_no2_x0_type0_noeri &
  (sa, ia, sa0, ia0, W2, d1, fc1, nir, nsym, psym, flops)


end subroutine g_if_sigma_ccov_ccoo_no2_x0_type0_noeri



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
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ccoo_no2_x0_type0_noeri &
  (s_a, i_a, s_a0, i_a0, W2_, D1_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
real(kind=8)                   :: W2_

! Intermediate arrays                
real*8, allocatable :: Z1_(:)
real*8, allocatable :: Z2_(:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1
! W2(a0,a) += (    1.00000000) D1(a1,a0) Fc1(a1,a) 
do s_a1 = 0, nir-1
if( &
IEOR(s_a0,s_a) == 0 .and. & 
IEOR(s_a1,s_a0) == 0 .and. &
IEOR(s_a1,s_a) == 0) then

if(psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D1(a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do
! Z2 <-- Fc1(a1,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1) =  &
  Fc1_(s_a, s_a1)%array(i_a, i_a1)
end do

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', 1,&
                     1,&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     1,&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     1)

! W2(a0,a)  <-- Z3
W2_ = &
    W2_ &
  + Z3_

! Flop count
flops = flops + 1 * &
                1 * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no2_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no2_x1_type0_noeri &
  (sa, ia, sa0, ia0, T2, W2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), W2, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no2_x1_type0_noeri &
  (sa, ia, sa0, ia0, av2_i, W2, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_ccov_ccoo_no2_x1_type0_noeri



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
subroutine g_sigma_ccov_ccoo_no2_x1_type0_noeri &
  (s_a, i_a, s_a0, i_a0, T2_, W2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W2_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_i, i_i
! S2(w,x,i,a) += (    2.00000000) T2(x,w,i,a0) W2(a0,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_i,s_a0) .and. &
IEOR(s_a0,s_a) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0) then

! Z1 <-- T2(x,w,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_w, i_i) =  &
  T2_(s_i, s_w, s_x)%array(i_i, i_w, i_x)
end do
end do
end do
! Z2 <-- W2(a0,a) 
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
                     2.00000000d+00, &
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

end subroutine g_sigma_ccov_ccoo_no2_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no3_x0_type0_noeri &
  (sa1, ia1, T2, W3, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa1, ia1
real(kind=8), intent(inout) :: T2(*), W3(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa1

call set_symblock_Xcca(sleft, W3, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_ccov_ccoo_no3_x0_type0_noeri &
  (sa1, ia1, av2_i, Xcca, d1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no3_x0_type0_noeri



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
subroutine g_sigma_ccov_ccoo_no3_x0_type0_noeri &
  (s_a1, i_a1, T2_, W3_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W3_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a0, i_a0, s_i, i_i
! W3(x,w,i,a1) += (    1.00000000) T2(x,w,a0,a1) D1(i,a0) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_x,s_w) == IEOR(s_i,s_a1) .and. & 
IEOR(s_x,s_w) == IEOR(s_a0,s_a1) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
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
! Z2 <-- T2(x,w,a0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x, i_w) =  &
  T2_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x)
end do
end do
end do

! Z3 <-- W3(x,w,i,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! W3(x,w,i,a1)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W3_(s_i, s_w, s_x)%array(i_i, i_w, i_x) = &
    W3_(s_i, s_w, s_x)%array(i_i, i_w, i_x) &
  + Z3_(i_i, i_x, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no3_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no3_x1_type0_noeri &
  (sa, ia, sa1, ia1, W3, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: W3(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa1

call set_symblock_Xcca(sleft, W3, nir, nsym, psym) ! -> Xcca (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no3_x1_type0_noeri &
  (sa, ia, sa1, ia1, Xcca, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no3_x1_type0_noeri



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
subroutine g_sigma_ccov_ccoo_no3_x1_type0_noeri &
  (s_a, i_a, s_a1, i_a1, W3_, S2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W3_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_i, i_i
! S2(w,x,i,a) += (    2.00000000) Fc1(a1,a) W3(x,w,i,a1) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a1,s_a) == 0 .and. &
IEOR(s_x,s_w) == IEOR(s_i,s_a1)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0) then

! Z1 <-- W3(x,w,i,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_w, i_i) =  &
  W3_(s_i, s_w, s_x)%array(i_i, i_w, i_x)
end do
end do
end do
! Z2 <-- Fc1(a1,a) 
Z2_ =  &
  Fc1_(s_a, s_a1)%array(i_a, i_a1)

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     1,&
                     2.00000000d+00, &
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

end subroutine g_sigma_ccov_ccoo_no3_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no4_x0_type0_noeri &
  (sa0, ia0, T2, W4, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0
real(kind=8), intent(inout) :: T2(*), W4(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xccaa(sleft, W4, nir, nsym, psym) ! -> Xccaa (allocate) 
call g_sigma_ccov_ccoo_no4_x0_type0_noeri &
  (sa0, ia0, av2_i, Xccaa, d1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xccaa)

end subroutine g_if_sigma_ccov_ccoo_no4_x0_type0_noeri



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
subroutine g_sigma_ccov_ccoo_no4_x0_type0_noeri &
  (s_a0, i_a0, T2_, W4_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a1, i_a1, s_i, i_i
! W4(x,w,i,a1) += (    1.00000000) T2(x,w,a1,a0) D1(i,a0) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_x,s_w) == IEOR(s_i,s_a1) .and. & 
IEOR(s_x,s_w) == IEOR(s_a1,s_a0) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0) then

! Z1 <-- T2(x,w,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_w, i_a1) =  &
  T2_(s_a1, s_w, s_x)%array(i_a1, i_w, i_x)
end do
end do
end do
! Z2 <-- D1(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z2_(i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do

! Z3 <-- W4(x,w,i,a1) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_O, s_i),&
                     1,&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1))

! W4(x,w,i,a1)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W4_(s_a1, s_i, s_w, s_x)%array(i_a1, i_i, i_w, i_x) = &
    W4_(s_a1, s_i, s_w, s_x)%array(i_a1, i_i, i_w, i_x) &
  + Z3_(i_x, i_w, i_a1, i_i)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_O, s_i) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no4_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no4_x1_type0_noeri &
  (sa, ia, W4, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: W4(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xccaa(sleft, W4, nir, nsym, psym) ! -> Xccaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no4_x1_type0_noeri &
  (sa, ia, Xccaa, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xccaa)

end subroutine g_if_sigma_ccov_ccoo_no4_x1_type0_noeri



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
subroutine g_sigma_ccov_ccoo_no4_x1_type0_noeri &
  (s_a, i_a, W4_, S2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_x, i_x, s_w, i_w, s_i, i_i
! S2(w,x,i,a) += (   -4.00000000) Fc1(a1,a) W4(x,w,i,a1) 
do s_a1 = 0, nir-1
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a1,s_a) == 0 .and. &
IEOR(s_x,s_w) == IEOR(s_i,s_a1)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- W4(x,w,i,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_w, i_i, i_a1) =  &
  W4_(s_a1, s_i, s_w, s_x)%array(i_a1, i_i, i_w, i_x)
end do
end do
end do
end do
! Z2 <-- Fc1(a1,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1) =  &
  Fc1_(s_a, s_a1)%array(i_a, i_a1)
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1),&
                     - 4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
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
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no4_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no5_x0_type0_noeri &
  (sa, ia, si, ii, T2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, si, ii
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(si, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no5_x0_type0_noeri &
  (sa, ia, si, ii, av2_i, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_ccov_ccoo_no5_x0_type0_noeri



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
subroutine g_sigma_ccov_ccoo_no5_x0_type0_noeri &
  (s_a, i_a, s_i, i_i, T2_, S2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a0, i_a0
! S2(w,x,i,a) += (    8.00000000) T2(x,w,a0,i) Fc1(a0,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_a0,s_i) .and. &
IEOR(s_a0,s_a) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(x,w,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_w, i_a0) =  &
  T2_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x)
end do
end do
end do
! Z2 <-- Fc1(a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0) =  &
  Fc1_(s_a, s_a0)%array(i_a, i_a0)
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w),&
                     1,&
                     psym(I_LENGTH,I_O, s_a0),&
                     8.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) * &
                1 * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no5_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no6_x0_type0_noeri &
  (sa, ia, sa0, ia0, T2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no6_x0_type0_noeri &
  (sa, ia, sa0, ia0, av2_i, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_ccov_ccoo_no6_x0_type0_noeri



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
subroutine g_sigma_ccov_ccoo_no6_x0_type0_noeri &
  (s_a, i_a, s_a0, i_a0, T2_, S2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_i, i_i
! S2(w,x,i,a) += (   -4.00000000) T2(x,w,i,a0) Fc1(a0,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_i,s_a0) .and. &
IEOR(s_a0,s_a) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0) then

! Z1 <-- T2(x,w,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_w, i_i) =  &
  T2_(s_i, s_w, s_x)%array(i_i, i_w, i_x)
end do
end do
end do
! Z2 <-- Fc1(a0,a) 
Z2_ =  &
  Fc1_(s_a, s_a0)%array(i_a, i_a0)

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     1,&
                     - 4.00000000d+00, &
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

end subroutine g_sigma_ccov_ccoo_no6_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no0_x0_type1_eri_o &
  (sa, ia, sa1, ia1, V2, W25, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: V2(*), W25(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa1,sa)

call set_symblock_Xaa(sleft, W25, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_ccov_ccoo_no0_x0_type1_eri_o &
  (sa, ia, sa1, ia1, h2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no0_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no0_x0_type1_eri_o &
  (s_a, i_a, s_a1, i_a1, V2_, W25_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
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

integer :: s_a3, i_a3, s_a2, i_a2, s_i, i_i, s_a0, i_a0
! W25(i,a0,a1,a) += (    1.00000000) V2(a1,a3,a2,a) D2(i,a0,a2,a3) 
do s_a3 = 0, nir-1
do s_a2 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_a1,s_a) .and. & 
IEOR(s_a1,s_a3) == IEOR(s_a2,s_a) .and. &
IEOR(s_i,s_a0) == IEOR(s_a2,s_a3)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3) > 0) then

! Z1 <-- D2(i,a0,a2,a3) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a2, i_a3) =  &
  D2_(s_a3, s_a2, s_a0, s_i)%array(i_a3, i_a2, i_a0, i_i)
end do
end do
end do
end do
! Z2 <-- V2(a1,a3,a2,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a3) =  &
  V2_(s_a, s_a2, s_a3)%array(i_a, i_a2, i_a3)
end do
end do

! Z3 <-- W25(i,a0,a1,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W25(i,a0,a1,a)  <-- Z3
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
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no0_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no0_x1_type1_eri_o &
  (sa, ia, sa1, ia1, T2, W25, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), W25(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa1,sa)

call set_symblock_Xaa(sleft, W25, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no0_x1_type1_eri_o &
  (sa, ia, sa1, ia1, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no0_x1_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no0_x1_type1_eri_o &
  (s_a, i_a, s_a1, i_a1, T2_, W25_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W25_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_x, i_x, s_a0, i_a0, s_i, i_i
! S2(w,x,i,a) += (    2.00000000) T2(w,x,a0,a1) W25(i,a0,a1,a) 
do s_w = 0, nir-1
do s_x = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_w,s_x) == IEOR(s_a0,s_a1) .and. &
IEOR(s_i,s_a0) == IEOR(s_a1,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W25(i,a0,a1,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0) =  &
  W25_(s_a0, s_i)%array(i_a0, i_i)
end do
end do
! Z2 <-- T2(w,x,a0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w, i_x) =  &
  T2_(s_a0, s_x, s_w)%array(i_a0, i_x, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_i, i_w, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no0_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no1_x0_type1_eri_o &
  (sa, ia, sa1, ia1, V2, W26, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: V2(*), W26(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa1,sa)

call set_symblock_Xaa(sleft, W26, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_ccov_ccoo_no1_x0_type1_eri_o &
  (sa, ia, sa1, ia1, h2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no1_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no1_x0_type1_eri_o &
  (s_a, i_a, s_a1, i_a1, V2_, W26_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W26_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_a2, i_a2, s_i, i_i, s_a0, i_a0
! W26(i,a0,a1,a) += (    1.00000000) V2(a1,a3,a2,a) D2(i,a3,a2,a0) 
do s_a3 = 0, nir-1
do s_a2 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_a1,s_a) .and. & 
IEOR(s_a1,s_a3) == IEOR(s_a2,s_a) .and. &
IEOR(s_i,s_a3) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- D2(i,a3,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a3, i_a2) =  &
  D2_(s_a0, s_a2, s_a3, s_i)%array(i_a0, i_a2, i_a3, i_i)
end do
end do
end do
end do
! Z2 <-- V2(a1,a3,a2,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a2) =  &
  V2_(s_a, s_a2, s_a3)%array(i_a, i_a2, i_a3)
end do
end do

! Z3 <-- W26(i,a0,a1,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W26(i,a0,a1,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W26_(s_a0, s_i)%array(i_a0, i_i) = &
    W26_(s_a0, s_i)%array(i_a0, i_i) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no1_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no1_x1_type1_eri_o &
  (sa, ia, sa1, ia1, T2, W26, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), W26(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa1,sa)

call set_symblock_Xaa(sleft, W26, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no1_x1_type1_eri_o &
  (sa, ia, sa1, ia1, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no1_x1_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no1_x1_type1_eri_o &
  (s_a, i_a, s_a1, i_a1, T2_, W26_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W26_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a0, i_a0, s_i, i_i
! S2(w,x,i,a) += (    2.00000000) T2(x,w,a0,a1) W26(i,a0,a1,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_a0,s_a1) .and. &
IEOR(s_i,s_a0) == IEOR(s_a1,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W26(i,a0,a1,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0) =  &
  W26_(s_a0, s_i)%array(i_a0, i_i)
end do
end do
! Z2 <-- T2(x,w,a0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x, i_w) =  &
  T2_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x)
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
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
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
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no1_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no2_x0_type1_eri_o &
  (sa, ia, sa2, ia2, V2, W31, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa2, ia2
real(kind=8), intent(inout) :: V2(*), W31(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa2, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa2,sa)

call set_symblock_Xaa(sleft, W31, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_ccov_ccoo_no2_x0_type1_eri_o &
  (sa, ia, sa2, ia2, h2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no2_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no2_x0_type1_eri_o &
  (s_a, i_a, s_a2, i_a2, V2_, W31_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a2, s_a2
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W31_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W31(i,a0,a2,a) += (    1.00000000) V2(a2,a,a3,a1) D2(i,a0,a3,a1) 
do s_a3 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_a2,s_a) .and. & 
IEOR(s_a2,s_a) == IEOR(s_a3,s_a1) .and. &
IEOR(s_i,s_a0) == IEOR(s_a3,s_a1)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(i,a0,a3,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a3, i_a1) =  &
  D2_(s_a1, s_a3, s_a0, s_i)%array(i_a1, i_a3, i_a0, i_i)
end do
end do
end do
end do
! Z2 <-- V2(a2,a,a3,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a1) =  &
  V2_(s_a1, s_a3, s_a)%array(i_a1, i_a3, i_a)
end do
end do

! Z3 <-- W31(i,a0,a2,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W31(i,a0,a2,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W31_(s_a0, s_i)%array(i_a0, i_i) = &
    W31_(s_a0, s_i)%array(i_a0, i_i) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no2_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no2_x1_type1_eri_o &
  (sa, ia, sa2, ia2, T2, W31, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa2, ia2
real(kind=8), intent(inout) :: T2(*), W31(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa2, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa2,sa)

call set_symblock_Xaa(sleft, W31, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no2_x1_type1_eri_o &
  (sa, ia, sa2, ia2, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no2_x1_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no2_x1_type1_eri_o &
  (s_a, i_a, s_a2, i_a2, T2_, W31_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a2, s_a2
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W31_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_x, i_x, s_a0, i_a0, s_i, i_i
! S2(w,x,i,a) += (   -4.00000000) T2(w,x,a0,a2) W31(i,a0,a2,a) 
do s_w = 0, nir-1
do s_x = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_w,s_x) == IEOR(s_a0,s_a2) .and. &
IEOR(s_i,s_a0) == IEOR(s_a2,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W31(i,a0,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0) =  &
  W31_(s_a0, s_i)%array(i_a0, i_i)
end do
end do
! Z2 <-- T2(w,x,a0,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w, i_x) =  &
  T2_(s_a0, s_x, s_w)%array(i_a0, i_x, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_i, i_w, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no2_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no3_x0_type1_eri_o &
  (sa, ia, sa2, ia2, V2, W32, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa2, ia2
real(kind=8), intent(inout) :: V2(*), W32(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa2, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa2,sa)

call set_symblock_Xaa(sleft, W32, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_ccov_ccoo_no3_x0_type1_eri_o &
  (sa, ia, sa2, ia2, h2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no3_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no3_x0_type1_eri_o &
  (s_a, i_a, s_a2, i_a2, V2_, W32_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a2, s_a2
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W32_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W32(i,a0,a2,a) += (    1.00000000) V2(a2,a,a3,a1) D2(i,a0,a3,a1) 
do s_a3 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_a2,s_a) .and. & 
IEOR(s_a2,s_a) == IEOR(s_a3,s_a1) .and. &
IEOR(s_i,s_a0) == IEOR(s_a3,s_a1)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(i,a0,a3,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a3, i_a1) =  &
  D2_(s_a1, s_a3, s_a0, s_i)%array(i_a1, i_a3, i_a0, i_i)
end do
end do
end do
end do
! Z2 <-- V2(a2,a,a3,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a1) =  &
  V2_(s_a1, s_a3, s_a)%array(i_a1, i_a3, i_a)
end do
end do

! Z3 <-- W32(i,a0,a2,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W32(i,a0,a2,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W32_(s_a0, s_i)%array(i_a0, i_i) = &
    W32_(s_a0, s_i)%array(i_a0, i_i) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no3_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no3_x1_type1_eri_o &
  (sa, ia, sa2, ia2, T2, W32, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa2, ia2
real(kind=8), intent(inout) :: T2(*), W32(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa2, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa2,sa)

call set_symblock_Xaa(sleft, W32, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no3_x1_type1_eri_o &
  (sa, ia, sa2, ia2, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no3_x1_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no3_x1_type1_eri_o &
  (s_a, i_a, s_a2, i_a2, T2_, W32_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a2, s_a2
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W32_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a0, i_a0, s_i, i_i
! S2(w,x,i,a) += (    2.00000000) T2(x,w,a0,a2) W32(i,a0,a2,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_a0,s_a2) .and. &
IEOR(s_i,s_a0) == IEOR(s_a2,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W32(i,a0,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0) =  &
  W32_(s_a0, s_i)%array(i_a0, i_i)
end do
end do
! Z2 <-- T2(x,w,a0,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x, i_w) =  &
  T2_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x)
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
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
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
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no3_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no4_x0_type1_eri_o &
  (sa, ia, sa2, ia2, V2, W33, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa2, ia2
real(kind=8), intent(inout) :: V2(*), W33(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa2, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa2,sa)

call set_symblock_Xaa(sleft, W33, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_ccov_ccoo_no4_x0_type1_eri_o &
  (sa, ia, sa2, ia2, h2_i, Xaa, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no4_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no4_x0_type1_eri_o &
  (s_a, i_a, s_a2, i_a2, V2_, W33_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a2, s_a2
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W33_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a1, i_a1, s_a0, i_a0
! W33(a0,i,a2,a) += (    1.00000000) V2(a2,a,i,a1) D1(a1,a0) 
do s_i = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_a0,s_i) == IEOR(s_a2,s_a) .and. & 
IEOR(s_a2,s_a) == IEOR(s_i,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(a2,a,i,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a1) =  &
  V2_(s_a1, s_i, s_a)%array(i_a1, i_i, i_a)
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

! Z3 <-- W33(a0,i,a2,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! W33(a0,i,a2,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W33_(s_i, s_a0)%array(i_i, i_a0) = &
    W33_(s_i, s_a0)%array(i_i, i_a0) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no4_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no4_x1_type1_eri_o &
  (sa, ia, sa2, ia2, T2, W33, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa2, ia2
real(kind=8), intent(inout) :: T2(*), W33(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa2, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa2,sa)

call set_symblock_Xaa(sleft, W33, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no4_x1_type1_eri_o &
  (sa, ia, sa2, ia2, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no4_x1_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no4_x1_type1_eri_o &
  (s_a, i_a, s_a2, i_a2, T2_, W33_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a2, s_a2
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W33_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_x, i_x, s_a0, i_a0, s_i, i_i
! S2(w,x,i,a) += (   -4.00000000) T2(w,x,a0,a2) W33(a0,i,a2,a) 
do s_w = 0, nir-1
do s_x = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_w,s_x) == IEOR(s_a0,s_a2) .and. &
IEOR(s_a0,s_i) == IEOR(s_a2,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W33(a0,i,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0) =  &
  W33_(s_i, s_a0)%array(i_i, i_a0)
end do
end do
! Z2 <-- T2(w,x,a0,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w, i_x) =  &
  T2_(s_a0, s_x, s_w)%array(i_a0, i_x, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_i, i_w, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no4_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no5_x0_type1_eri_o &
  (sa, ia, sa2, ia2, V2, W34, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa2, ia2
real(kind=8), intent(inout) :: V2(*), W34(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa2, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa2,sa)

call set_symblock_Xaa(sleft, W34, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_ccov_ccoo_no5_x0_type1_eri_o &
  (sa, ia, sa2, ia2, h2_i, Xaa, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no5_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no5_x0_type1_eri_o &
  (s_a, i_a, s_a2, i_a2, V2_, W34_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a2, s_a2
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W34_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a1, i_a1, s_a0, i_a0
! W34(a0,i,a2,a) += (    1.00000000) V2(a2,a,i,a1) D1(a1,a0) 
do s_i = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_a0,s_i) == IEOR(s_a2,s_a) .and. & 
IEOR(s_a2,s_a) == IEOR(s_i,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(a2,a,i,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a1) =  &
  V2_(s_a1, s_i, s_a)%array(i_a1, i_i, i_a)
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

! Z3 <-- W34(a0,i,a2,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! W34(a0,i,a2,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W34_(s_i, s_a0)%array(i_i, i_a0) = &
    W34_(s_i, s_a0)%array(i_i, i_a0) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no5_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no5_x1_type1_eri_o &
  (sa, ia, sa2, ia2, T2, W34, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa2, ia2
real(kind=8), intent(inout) :: T2(*), W34(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa2, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa2,sa)

call set_symblock_Xaa(sleft, W34, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no5_x1_type1_eri_o &
  (sa, ia, sa2, ia2, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no5_x1_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no5_x1_type1_eri_o &
  (s_a, i_a, s_a2, i_a2, T2_, W34_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a2, s_a2
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W34_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a0, i_a0, s_i, i_i
! S2(w,x,i,a) += (    2.00000000) T2(x,w,a0,a2) W34(a0,i,a2,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_a0,s_a2) .and. &
IEOR(s_a0,s_i) == IEOR(s_a2,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W34(a0,i,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0) =  &
  W34_(s_i, s_a0)%array(i_i, i_a0)
end do
end do
! Z2 <-- T2(x,w,a0,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x, i_w) =  &
  T2_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x)
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
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
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
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no5_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no6_x0_type1_eri_o &
  (sa, ia, sa0, ia0, V2, W35, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: V2(*), W35
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa0,sa)

call g_sigma_ccov_ccoo_no6_x0_type1_eri_o &
  (sa, ia, sa0, ia0, h2_i, W35, d1, nir, nsym, psym, flops)

deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no6_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no6_x0_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, V2_, W35_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W35_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a1, i_a1
! W35(a0,a) += (    1.00000000) V2(a0,a2,a1,a) D1(a2,a1) 
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(s_a0,s_a) == 0 .and. & 
IEOR(s_a0,s_a2) == IEOR(s_a1,s_a) .and. &
IEOR(s_a2,s_a1) == 0) then

if(psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(a0,a2,a1,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_a1) =  &
  V2_(s_a, s_a1, s_a2)%array(i_a, i_a1, i_a2)
end do
end do
! Z2 <-- D1(a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1) =  &
  D1_(s_a1, s_a2)%array(i_a1, i_a2)
end do
end do

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', 1,&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     1,&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     1)

! W35(a0,a)  <-- Z3
W35_ = &
    W35_ &
  + Z3_

! Flop count
flops = flops + 1 * &
                1 * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no6_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no6_x1_type1_eri_o &
  (sa, ia, sa0, ia0, T2, W35, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), W35, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no6_x1_type1_eri_o &
  (sa, ia, sa0, ia0, av2_i, W35, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_ccov_ccoo_no6_x1_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no6_x1_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, T2_, W35_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W35_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_i, i_i
! S2(w,x,i,a) += (    2.00000000) T2(x,w,i,a0) W35(a0,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_i,s_a0) .and. &
IEOR(s_a0,s_a) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0) then

! Z1 <-- T2(x,w,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_w, i_i) =  &
  T2_(s_i, s_w, s_x)%array(i_i, i_w, i_x)
end do
end do
end do
! Z2 <-- W35(a0,a) 
Z2_ =  &
  W35_

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     1,&
                     2.00000000d+00, &
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

end subroutine g_sigma_ccov_ccoo_no6_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no7_x0_type1_eri_o &
  (sa, ia, sa0, ia0, V2, W36, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: V2(*), W36
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa0,sa)

call g_sigma_ccov_ccoo_no7_x0_type1_eri_o &
  (sa, ia, sa0, ia0, h2_i, W36, d1, nir, nsym, psym, flops)

deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no7_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no7_x0_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, V2_, W36_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W36_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a1, i_a1
! W36(a0,a) += (    1.00000000) V2(a0,a2,a1,a) D1(a2,a1) 
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(s_a0,s_a) == 0 .and. & 
IEOR(s_a0,s_a2) == IEOR(s_a1,s_a) .and. &
IEOR(s_a2,s_a1) == 0) then

if(psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(a0,a2,a1,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_a1) =  &
  V2_(s_a, s_a1, s_a2)%array(i_a, i_a1, i_a2)
end do
end do
! Z2 <-- D1(a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1) =  &
  D1_(s_a1, s_a2)%array(i_a1, i_a2)
end do
end do

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', 1,&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     1,&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     1)

! W36(a0,a)  <-- Z3
W36_ = &
    W36_ &
  + Z3_

! Flop count
flops = flops + 1 * &
                1 * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no7_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no7_x1_type1_eri_o &
  (sa, ia, sa0, ia0, T2, W36, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), W36, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no7_x1_type1_eri_o &
  (sa, ia, sa0, ia0, av2_i, W36, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_ccov_ccoo_no7_x1_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no7_x1_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, T2_, W36_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W36_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_x, i_x, s_i, i_i
! S2(w,x,i,a) += (   -4.00000000) T2(w,x,i,a0) W36(a0,a) 
do s_w = 0, nir-1
do s_x = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_w,s_x) == IEOR(s_i,s_a0) .and. &
IEOR(s_a0,s_a) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) > 0) then

! Z1 <-- T2(w,x,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_x, i_i) =  &
  T2_(s_i, s_x, s_w)%array(i_i, i_x, i_w)
end do
end do
end do
! Z2 <-- W36(a0,a) 
Z2_ =  &
  W36_

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     1,&
                     - 4.00000000d+00, &
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

end subroutine g_sigma_ccov_ccoo_no7_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no8_x0_type1_eri_o &
  (sa, ia, sa0, ia0, T2, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no8_x0_type1_eri_o &
  (sa, ia, sa0, ia0, av2_i, h2_i, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no8_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no8_x0_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_c0, i_c0, s_w, i_w, s_i, i_i
! S2(w,x,i,a) += (   -8.00000000) V2(a0,a,x,c0) T2(w,c0,i,a0) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a0,s_a) == IEOR(s_x,s_c0) .and. &
IEOR(s_w,s_c0) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(w,c0,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i, i_c0) =  &
  T2_(s_i, s_c0, s_w)%array(i_i, i_c0, i_w)
end do
end do
end do
! Z2 <-- V2(a0,a,x,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_x) =  &
  V2_(s_c0, s_x, s_a)%array(i_c0, i_x, i_a)
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_c0),&
                     - 8.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_w, i_i, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no8_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no9_x0_type1_eri_o &
  (sa, ia, sa0, ia0, T2, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no9_x0_type1_eri_o &
  (sa, ia, sa0, ia0, av2_i, h2_i, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no9_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no9_x0_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_x, i_x, s_i, i_i
! S2(w,x,i,a) += (    4.00000000) V2(a0,a,w,c0) T2(x,c0,i,a0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_x = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a0,s_a) == IEOR(s_w,s_c0) .and. &
IEOR(s_x,s_c0) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(x,c0,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_i, i_c0) =  &
  T2_(s_i, s_c0, s_x)%array(i_i, i_c0, i_x)
end do
end do
end do
! Z2 <-- V2(a0,a,w,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  V2_(s_c0, s_w, s_a)%array(i_c0, i_w, i_a)
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_ccov_ccoo_no9_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no10_x0_type1_eri_o &
  (sa, ia, sa0, ia0, T2, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no10_x0_type1_eri_o &
  (sa, ia, sa0, ia0, av2_i, h2_i, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no10_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no10_x0_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_c0, i_c0, s_w, i_w, s_i, i_i
! S2(w,x,i,a) += (    4.00000000) V2(a0,a,x,c0) T2(c0,w,i,a0) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a0,s_a) == IEOR(s_x,s_c0) .and. &
IEOR(s_c0,s_w) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c0,w,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i, i_c0) =  &
  T2_(s_i, s_w, s_c0)%array(i_i, i_w, i_c0)
end do
end do
end do
! Z2 <-- V2(a0,a,x,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_x) =  &
  V2_(s_c0, s_x, s_a)%array(i_c0, i_x, i_a)
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_c0),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_w, i_i, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no10_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no11_x0_type1_eri_o &
  (sa, ia, sa0, ia0, T2, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no11_x0_type1_eri_o &
  (sa, ia, sa0, ia0, av2_i, h2_i, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no11_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no11_x0_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_x, i_x, s_i, i_i
! S2(w,x,i,a) += (   -8.00000000) V2(a0,a,w,c0) T2(c0,x,i,a0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_x = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a0,s_a) == IEOR(s_w,s_c0) .and. &
IEOR(s_c0,s_x) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c0,x,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_i, i_c0) =  &
  T2_(s_i, s_x, s_c0)%array(i_i, i_x, i_c0)
end do
end do
end do
! Z2 <-- V2(a0,a,w,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  V2_(s_c0, s_w, s_a)%array(i_c0, i_w, i_a)
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     - 8.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_ccov_ccoo_no11_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no12_x0_type1_eri_o &
  (sa, ia, sa1, ia1, V2, W37, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: V2(*), W37
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa1,sa)

call g_sigma_ccov_ccoo_no12_x0_type1_eri_o &
  (sa, ia, sa1, ia1, h2_i, W37, d1, nir, nsym, psym, flops)

deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no12_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no12_x0_type1_eri_o &
  (s_a, i_a, s_a1, i_a1, V2_, W37_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W37_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a0, i_a0
! W37(a1,a) += (    1.00000000) V2(a1,a,a2,a0) D1(a2,a0) 
do s_a2 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_a1,s_a) == 0 .and. & 
IEOR(s_a1,s_a) == IEOR(s_a2,s_a0) .and. &
IEOR(s_a2,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(a1,a,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_a0) =  &
  V2_(s_a0, s_a2, s_a)%array(i_a0, i_a2, i_a)
end do
end do
! Z2 <-- D1(a2,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a0) =  &
  D1_(s_a0, s_a2)%array(i_a0, i_a2)
end do
end do

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', 1,&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     1,&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     1)

! W37(a1,a)  <-- Z3
W37_ = &
    W37_ &
  + Z3_

! Flop count
flops = flops + 1 * &
                1 * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no12_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no12_x1_type1_eri_o &
  (sa, ia, sa1, ia1, T2, W37, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), W37, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa1,sa)

call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no12_x1_type1_eri_o &
  (sa, ia, sa1, ia1, av2_i, W37, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_ccov_ccoo_no12_x1_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no12_x1_type1_eri_o &
  (s_a, i_a, s_a1, i_a1, T2_, W37_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W37_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_x, i_x, s_i, i_i
! S2(w,x,i,a) += (    8.00000000) T2(w,x,i,a1) W37(a1,a) 
do s_w = 0, nir-1
do s_x = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_w,s_x) == IEOR(s_i,s_a1) .and. &
IEOR(s_a1,s_a) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) > 0) then

! Z1 <-- T2(w,x,i,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_x, i_i) =  &
  T2_(s_i, s_x, s_w)%array(i_i, i_x, i_w)
end do
end do
end do
! Z2 <-- W37(a1,a) 
Z2_ =  &
  W37_

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     1,&
                     8.00000000d+00, &
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

end subroutine g_sigma_ccov_ccoo_no12_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no13_x0_type1_eri_o &
  (sa, ia, sa1, ia1, V2, W38, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: V2(*), W38
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa1,sa)

call g_sigma_ccov_ccoo_no13_x0_type1_eri_o &
  (sa, ia, sa1, ia1, h2_i, W38, d1, nir, nsym, psym, flops)

deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no13_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no13_x0_type1_eri_o &
  (s_a, i_a, s_a1, i_a1, V2_, W38_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W38_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a0, i_a0
! W38(a1,a) += (    1.00000000) V2(a1,a,a2,a0) D1(a2,a0) 
do s_a2 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_a1,s_a) == 0 .and. & 
IEOR(s_a1,s_a) == IEOR(s_a2,s_a0) .and. &
IEOR(s_a2,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(a1,a,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_a0) =  &
  V2_(s_a0, s_a2, s_a)%array(i_a0, i_a2, i_a)
end do
end do
! Z2 <-- D1(a2,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a0) =  &
  D1_(s_a0, s_a2)%array(i_a0, i_a2)
end do
end do

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', 1,&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     1,&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     1)

! W38(a1,a)  <-- Z3
W38_ = &
    W38_ &
  + Z3_

! Flop count
flops = flops + 1 * &
                1 * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no13_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no13_x1_type1_eri_o &
  (sa, ia, sa1, ia1, T2, W38, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), W38, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa1,sa)

call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no13_x1_type1_eri_o &
  (sa, ia, sa1, ia1, av2_i, W38, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_ccov_ccoo_no13_x1_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no13_x1_type1_eri_o &
  (s_a, i_a, s_a1, i_a1, T2_, W38_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W38_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_i, i_i
! S2(w,x,i,a) += (   -4.00000000) T2(x,w,i,a1) W38(a1,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_i,s_a1) .and. &
IEOR(s_a1,s_a) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0) then

! Z1 <-- T2(x,w,i,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_w, i_i) =  &
  T2_(s_i, s_w, s_x)%array(i_i, i_w, i_x)
end do
end do
end do
! Z2 <-- W38(a1,a) 
Z2_ =  &
  W38_

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     1,&
                     - 4.00000000d+00, &
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

end subroutine g_sigma_ccov_ccoo_no13_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no14_x0_type1_eri_o &
  (sa, ia, sa1, ia1, V2, W39, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: V2(*), W39(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa1,sa)

call set_symblock_Xaa(sleft, W39, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_ccov_ccoo_no14_x0_type1_eri_o &
  (sa, ia, sa1, ia1, h2_i, Xaa, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no14_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no14_x0_type1_eri_o &
  (s_a, i_a, s_a1, i_a1, V2_, W39_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W39_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a2, i_a2, s_a0, i_a0
! W39(a0,i,a1,a) += (    1.00000000) V2(a1,i,a2,a) D1(a2,a0) 
do s_i = 0, nir-1
do s_a2 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_a0,s_i) == IEOR(s_a1,s_a) .and. & 
IEOR(s_a1,s_i) == IEOR(s_a2,s_a) .and. &
IEOR(s_a2,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(a1,i,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a2) =  &
  V2_(s_a, s_a2, s_i)%array(i_a, i_a2, i_i)
end do
end do
! Z2 <-- D1(a2,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a0) =  &
  D1_(s_a0, s_a2)%array(i_a0, i_a2)
end do
end do

! Z3 <-- W39(a0,i,a1,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! W39(a0,i,a1,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W39_(s_i, s_a0)%array(i_i, i_a0) = &
    W39_(s_i, s_a0)%array(i_i, i_a0) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no14_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no14_x1_type1_eri_o &
  (sa, ia, sa1, ia1, T2, W39, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), W39(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa1,sa)

call set_symblock_Xaa(sleft, W39, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no14_x1_type1_eri_o &
  (sa, ia, sa1, ia1, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no14_x1_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no14_x1_type1_eri_o &
  (s_a, i_a, s_a1, i_a1, T2_, W39_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W39_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_x, i_x, s_a0, i_a0, s_i, i_i
! S2(w,x,i,a) += (    2.00000000) T2(w,x,a0,a1) W39(a0,i,a1,a) 
do s_w = 0, nir-1
do s_x = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_w,s_x) == IEOR(s_a0,s_a1) .and. &
IEOR(s_a0,s_i) == IEOR(s_a1,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W39(a0,i,a1,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0) =  &
  W39_(s_i, s_a0)%array(i_i, i_a0)
end do
end do
! Z2 <-- T2(w,x,a0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w, i_x) =  &
  T2_(s_a0, s_x, s_w)%array(i_a0, i_x, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_i, i_w, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no14_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no15_x0_type1_eri_o &
  (sa, ia, sa1, ia1, V2, W40, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: V2(*), W40(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa1,sa)

call set_symblock_Xaa(sleft, W40, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_ccov_ccoo_no15_x0_type1_eri_o &
  (sa, ia, sa1, ia1, h2_i, Xaa, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no15_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no15_x0_type1_eri_o &
  (s_a, i_a, s_a1, i_a1, V2_, W40_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W40_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a2, i_a2, s_a0, i_a0
! W40(a0,i,a1,a) += (    1.00000000) V2(a1,i,a2,a) D1(a2,a0) 
do s_i = 0, nir-1
do s_a2 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_a0,s_i) == IEOR(s_a1,s_a) .and. & 
IEOR(s_a1,s_i) == IEOR(s_a2,s_a) .and. &
IEOR(s_a2,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(a1,i,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a2) =  &
  V2_(s_a, s_a2, s_i)%array(i_a, i_a2, i_i)
end do
end do
! Z2 <-- D1(a2,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a0) =  &
  D1_(s_a0, s_a2)%array(i_a0, i_a2)
end do
end do

! Z3 <-- W40(a0,i,a1,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! W40(a0,i,a1,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W40_(s_i, s_a0)%array(i_i, i_a0) = &
    W40_(s_i, s_a0)%array(i_i, i_a0) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no15_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no15_x1_type1_eri_o &
  (sa, ia, sa1, ia1, T2, W40, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), W40(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa1,sa)

call set_symblock_Xaa(sleft, W40, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no15_x1_type1_eri_o &
  (sa, ia, sa1, ia1, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no15_x1_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no15_x1_type1_eri_o &
  (s_a, i_a, s_a1, i_a1, T2_, W40_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W40_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a0, i_a0, s_i, i_i
! S2(w,x,i,a) += (   -4.00000000) T2(x,w,a0,a1) W40(a0,i,a1,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_a0,s_a1) .and. &
IEOR(s_a0,s_i) == IEOR(s_a1,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W40(a0,i,a1,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0) =  &
  W40_(s_i, s_a0)%array(i_i, i_a0)
end do
end do
! Z2 <-- T2(x,w,a0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x, i_w) =  &
  T2_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x)
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
                     psym(I_LENGTH,I_O, s_a0),&
                     - 4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
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
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no15_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no16_x0_type1_eri_o &
  (sa, ia, sa0, ia0, V2, W41, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: V2(*), W41(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xaa(sleft, W41, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_ccov_ccoo_no16_x0_type1_eri_o &
  (sa, ia, sa0, ia0, h2_i, Xaa, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no16_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no16_x0_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, V2_, W41_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W41_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a1, i_a1, s_i, i_i
! W41(i,a1,a0,a) += (    1.00000000) V2(a0,a2,a1,a) D1(i,a2) 
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_a1) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a0,s_a2) == IEOR(s_a1,s_a) .and. &
IEOR(s_i,s_a2) == 0) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(a0,a2,a1,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_a2) =  &
  V2_(s_a, s_a1, s_a2)%array(i_a, i_a1, i_a2)
end do
end do
! Z2 <-- D1(i,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_i) =  &
  D1_(s_a2, s_i)%array(i_a2, i_i)
end do
end do

! Z3 <-- W41(i,a1,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W41(i,a1,a0,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W41_(s_a1, s_i)%array(i_a1, i_i) = &
    W41_(s_a1, s_i)%array(i_a1, i_i) &
  + Z3_(i_a1, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no16_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no16_x1_type1_eri_o &
  (sa, ia, sa0, ia0, T2, W41, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), W41(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xaa(sleft, W41, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no16_x1_type1_eri_o &
  (sa, ia, sa0, ia0, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no16_x1_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no16_x1_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, T2_, W41_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W41_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a1, i_a1, s_i, i_i
! S2(w,x,i,a) += (   -4.00000000) T2(x,w,a1,a0) W41(i,a1,a0,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_a1,s_a0) .and. &
IEOR(s_i,s_a1) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- W41(i,a1,a0,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a1) =  &
  W41_(s_a1, s_i)%array(i_a1, i_i)
end do
end do
! Z2 <-- T2(x,w,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_x, i_w) =  &
  T2_(s_a1, s_w, s_x)%array(i_a1, i_w, i_x)
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
                     psym(I_LENGTH,I_O, s_a1),&
                     - 4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
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
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no16_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no17_x0_type1_eri_o &
  (sa, ia, sa0, ia0, V2, W42, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: V2(*), W42(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xaa(sleft, W42, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_ccov_ccoo_no17_x0_type1_eri_o &
  (sa, ia, sa0, ia0, h2_i, Xaa, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no17_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no17_x0_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, V2_, W42_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W42_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a1, i_a1, s_i, i_i
! W42(i,a1,a0,a) += (    1.00000000) V2(a0,a2,a1,a) D1(i,a2) 
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_a1) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a0,s_a2) == IEOR(s_a1,s_a) .and. &
IEOR(s_i,s_a2) == 0) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(a0,a2,a1,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_a2) =  &
  V2_(s_a, s_a1, s_a2)%array(i_a, i_a1, i_a2)
end do
end do
! Z2 <-- D1(i,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_i) =  &
  D1_(s_a2, s_i)%array(i_a2, i_i)
end do
end do

! Z3 <-- W42(i,a1,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W42(i,a1,a0,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W42_(s_a1, s_i)%array(i_a1, i_i) = &
    W42_(s_a1, s_i)%array(i_a1, i_i) &
  + Z3_(i_a1, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no17_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no17_x1_type1_eri_o &
  (sa, ia, sa0, ia0, T2, W42, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), W42(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xaa(sleft, W42, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no17_x1_type1_eri_o &
  (sa, ia, sa0, ia0, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no17_x1_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no17_x1_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, T2_, W42_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W42_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_x, i_x, s_a1, i_a1, s_i, i_i
! S2(w,x,i,a) += (    2.00000000) T2(w,x,a1,a0) W42(i,a1,a0,a) 
do s_w = 0, nir-1
do s_x = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_w,s_x) == IEOR(s_a1,s_a0) .and. &
IEOR(s_i,s_a1) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- W42(i,a1,a0,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a1) =  &
  W42_(s_a1, s_i)%array(i_a1, i_i)
end do
end do
! Z2 <-- T2(w,x,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w, i_x) =  &
  T2_(s_a1, s_x, s_w)%array(i_a1, i_x, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a1),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_i, i_w, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no17_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no18_x0_type1_eri_o &
  (sa, ia, sa0, ia0, T2, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no18_x0_type1_eri_o &
  (sa, ia, sa0, ia0, av2_i, h2_i, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no18_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no18_x0_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a1, i_a1, s_x, i_x, s_w, i_w
! S2(w,x,i,a) += (    8.00000000) V2(a0,i,a1,a) T2(x,w,a1,a0) 
do s_i = 0, nir-1
do s_a1 = 0, nir-1
do s_x = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a0,s_i) == IEOR(s_a1,s_a) .and. &
IEOR(s_x,s_w) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(a0,i,a1,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a1) =  &
  V2_(s_a, s_a1, s_i)%array(i_a, i_a1, i_i)
end do
end do
! Z2 <-- T2(x,w,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_x, i_w) =  &
  T2_(s_a1, s_w, s_x)%array(i_a1, i_w, i_x)
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
                     psym(I_LENGTH,I_O, s_a1),&
                     8.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
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
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no18_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no19_x0_type1_eri_o &
  (sa, ia, sa1, ia1, T2, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no19_x0_type1_eri_o &
  (sa, ia, sa1, ia1, av2_i, h2_i, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no19_x0_type1_eri_o



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
subroutine g_sigma_ccov_ccoo_no19_x0_type1_eri_o &
  (s_a, i_a, s_a1, i_a1, T2_, V2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_w, i_w, s_x, i_x
! S2(w,x,i,a) += (   -4.00000000) V2(a1,i,a0,a) T2(w,x,a0,a1) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a1,s_i) == IEOR(s_a0,s_a) .and. &
IEOR(s_w,s_x) == IEOR(s_a0,s_a1)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(a1,i,a0,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0) =  &
  V2_(s_a, s_a0, s_i)%array(i_a, i_a0, i_i)
end do
end do
! Z2 <-- T2(w,x,a0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w, i_x) =  &
  T2_(s_a0, s_x, s_w)%array(i_a0, i_x, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_i, i_w, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no19_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no0_x0_type0_eri_v &
  (sa0, ia0, T2, W5, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0
real(kind=8), intent(inout) :: T2(*), W5(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xccaa(sleft, W5, nir, nsym, psym) ! -> Xccaa (allocate) 
call g_sigma_ccov_ccoo_no0_x0_type0_eri_v &
  (sa0, ia0, av2_i, Xccaa, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xccaa)

end subroutine g_if_sigma_ccov_ccoo_no0_x0_type0_eri_v



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
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ccoo_no0_x0_type0_eri_v &
  (s_a0, i_a0, T2_, W5_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W5_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_a1, i_a1, s_i, i_i, s_a2, i_a2
! W5(w,c0,i,a2) += (    1.00000000) T2(w,c0,a1,a0) D2(i,a1,a2,a0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a2 = 0, nir-1
if( &
IEOR(s_w,s_c0) == IEOR(s_i,s_a2) .and. & 
IEOR(s_w,s_c0) == IEOR(s_a1,s_a0) .and. &
IEOR(s_i,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(i,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a2, i_a1) =  &
  D2_(s_a0, s_a2, s_a1, s_i)%array(i_a0, i_a2, i_a1, i_i)
end do
end do
end do
! Z2 <-- T2(w,c0,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w, i_c0) =  &
  T2_(s_a1, s_c0, s_w)%array(i_a1, i_c0, i_w)
end do
end do
end do

! Z3 <-- W5(w,c0,i,a2) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2))

! W5(w,c0,i,a2)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W5_(s_a2, s_i, s_c0, s_w)%array(i_a2, i_i, i_c0, i_w) = &
    W5_(s_a2, s_i, s_c0, s_w)%array(i_a2, i_i, i_c0, i_w) &
  + Z3_(i_i, i_a2, i_w, i_c0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no0_x0_type0_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no1_x0_type0_eri_v &
  (sa0, ia0, T2, W6, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0
real(kind=8), intent(inout) :: T2(*), W6(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xccaa(sleft, W6, nir, nsym, psym) ! -> Xccaa (allocate) 
call g_sigma_ccov_ccoo_no1_x0_type0_eri_v &
  (sa0, ia0, av2_i, Xccaa, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xccaa)

end subroutine g_if_sigma_ccov_ccoo_no1_x0_type0_eri_v



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
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ccoo_no1_x0_type0_eri_v &
  (s_a0, i_a0, T2_, W6_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W6_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_c0, i_c0, s_a1, i_a1, s_i, i_i, s_a2, i_a2
! W6(x,c0,i,a2) += (    1.00000000) T2(x,c0,a1,a0) D2(i,a1,a2,a0) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a2 = 0, nir-1
if( &
IEOR(s_x,s_c0) == IEOR(s_i,s_a2) .and. & 
IEOR(s_x,s_c0) == IEOR(s_a1,s_a0) .and. &
IEOR(s_i,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(i,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a2, i_a1) =  &
  D2_(s_a0, s_a2, s_a1, s_i)%array(i_a0, i_a2, i_a1, i_i)
end do
end do
end do
! Z2 <-- T2(x,c0,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_x, i_c0) =  &
  T2_(s_a1, s_c0, s_x)%array(i_a1, i_c0, i_x)
end do
end do
end do

! Z3 <-- W6(x,c0,i,a2) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2))

! W6(x,c0,i,a2)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W6_(s_a2, s_i, s_c0, s_x)%array(i_a2, i_i, i_c0, i_x) = &
    W6_(s_a2, s_i, s_c0, s_x)%array(i_a2, i_i, i_c0, i_x) &
  + Z3_(i_i, i_a2, i_x, i_c0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) * &
                psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no1_x0_type0_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no2_x0_type0_eri_v &
  (sa0, ia0, T2, W7, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0
real(kind=8), intent(inout) :: T2(*), W7(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xccaa(sleft, W7, nir, nsym, psym) ! -> Xccaa (allocate) 
call g_sigma_ccov_ccoo_no2_x0_type0_eri_v &
  (sa0, ia0, av2_i, Xccaa, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xccaa)

end subroutine g_if_sigma_ccov_ccoo_no2_x0_type0_eri_v



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
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ccoo_no2_x0_type0_eri_v &
  (s_a0, i_a0, T2_, W7_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W7_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_a1, i_a1, s_i, i_i, s_a2, i_a2
! W7(w,c0,i,a2) += (    1.00000000) T2(w,c0,a1,a0) D2(i,a1,a2,a0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a2 = 0, nir-1
if( &
IEOR(s_w,s_c0) == IEOR(s_i,s_a2) .and. & 
IEOR(s_w,s_c0) == IEOR(s_a1,s_a0) .and. &
IEOR(s_i,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(i,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a2, i_a1) =  &
  D2_(s_a0, s_a2, s_a1, s_i)%array(i_a0, i_a2, i_a1, i_i)
end do
end do
end do
! Z2 <-- T2(w,c0,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w, i_c0) =  &
  T2_(s_a1, s_c0, s_w)%array(i_a1, i_c0, i_w)
end do
end do
end do

! Z3 <-- W7(w,c0,i,a2) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2))

! W7(w,c0,i,a2)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W7_(s_a2, s_i, s_c0, s_w)%array(i_a2, i_i, i_c0, i_w) = &
    W7_(s_a2, s_i, s_c0, s_w)%array(i_a2, i_i, i_c0, i_w) &
  + Z3_(i_i, i_a2, i_w, i_c0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no2_x0_type0_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no3_x0_type0_eri_v &
  (sa1, ia1, T2, W8, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa1, ia1
real(kind=8), intent(inout) :: T2(*), W8(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xccaa(sleft, W8, nir, nsym, psym) ! -> Xccaa (allocate) 
call g_sigma_ccov_ccoo_no3_x0_type0_eri_v &
  (sa1, ia1, av2_i, Xccaa, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xccaa)

end subroutine g_if_sigma_ccov_ccoo_no3_x0_type0_eri_v



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
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ccoo_no3_x0_type0_eri_v &
  (s_a1, i_a1, T2_, W8_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W8_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_c0, i_c0, s_a0, i_a0, s_i, i_i, s_a2, i_a2
! W8(x,c0,i,a2) += (    1.00000000) T2(x,c0,a0,a1) D2(i,a1,a2,a0) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_a2 = 0, nir-1
if( &
IEOR(s_x,s_c0) == IEOR(s_i,s_a2) .and. & 
IEOR(s_x,s_c0) == IEOR(s_a0,s_a1) .and. &
IEOR(s_i,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(i,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a2, i_a0) =  &
  D2_(s_a0, s_a2, s_a1, s_i)%array(i_a0, i_a2, i_a1, i_i)
end do
end do
end do
! Z2 <-- T2(x,c0,a0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x, i_c0) =  &
  T2_(s_a0, s_c0, s_x)%array(i_a0, i_c0, i_x)
end do
end do
end do

! Z3 <-- W8(x,c0,i,a2) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2))

! W8(x,c0,i,a2)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W8_(s_a2, s_i, s_c0, s_x)%array(i_a2, i_i, i_c0, i_x) = &
    W8_(s_a2, s_i, s_c0, s_x)%array(i_a2, i_i, i_c0, i_x) &
  + Z3_(i_i, i_a2, i_x, i_c0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) * &
                psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no3_x0_type0_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no0_x0_type1_eri_v &
  (sa, ia, V2, W5, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), W5(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_Xccaa(sleft, W5, nir, nsym, psym) ! -> Xccaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no0_x0_type1_eri_v &
  (sa, ia, h2_i, Xccaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)
deallocate(Xccaa)

end subroutine g_if_sigma_ccov_ccoo_no0_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no0_x0_type1_eri_v &
  (s_a, i_a, V2_, W5_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W5_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_c0, i_c0, s_a2, i_a2, s_w, i_w, s_i, i_i
! S2(w,x,i,a) += (    4.00000000) V2(a,x,c0,a2) W5(w,c0,i,a2) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_a2 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a,s_x) == IEOR(s_c0,s_a2) .and. &
IEOR(s_w,s_c0) == IEOR(s_i,s_a2)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- W5(w,c0,i,a2) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i, i_c0, i_a2) =  &
  W5_(s_a2, s_i, s_c0, s_w)%array(i_a2, i_i, i_c0, i_w)
end do
end do
end do
end do
! Z2 <-- V2(a,x,c0,a2) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a2, i_x) =  &
  V2_(s_a2, s_c0, s_x)%array(i_a2, i_c0, i_x)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_w, i_i, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no0_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no1_x0_type1_eri_v &
  (sa, ia, V2, W6, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), W6(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_Xccaa(sleft, W6, nir, nsym, psym) ! -> Xccaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no1_x0_type1_eri_v &
  (sa, ia, h2_i, Xccaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)
deallocate(Xccaa)

end subroutine g_if_sigma_ccov_ccoo_no1_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no1_x0_type1_eri_v &
  (s_a, i_a, V2_, W6_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W6_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_a2, i_a2, s_x, i_x, s_i, i_i
! S2(w,x,i,a) += (   -2.00000000) V2(a,w,c0,a2) W6(x,c0,i,a2) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a2 = 0, nir-1
do s_x = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_c0,s_a2) .and. &
IEOR(s_x,s_c0) == IEOR(s_i,s_a2)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- W6(x,c0,i,a2) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_i, i_c0, i_a2) =  &
  W6_(s_a2, s_i, s_c0, s_x)%array(i_a2, i_i, i_c0, i_x)
end do
end do
end do
end do
! Z2 <-- V2(a,w,c0,a2) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a2, i_w) =  &
  V2_(s_a2, s_c0, s_w)%array(i_a2, i_c0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no1_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no2_x0_type1_eri_v &
  (sa, ia, V2, W7, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), W7(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_Xccaa(sleft, W7, nir, nsym, psym) ! -> Xccaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no2_x0_type1_eri_v &
  (sa, ia, h2_i, Xccaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)
deallocate(Xccaa)

end subroutine g_if_sigma_ccov_ccoo_no2_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no2_x0_type1_eri_v &
  (s_a, i_a, V2_, W7_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W7_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_x, i_x, s_c0, i_c0, s_w, i_w, s_i, i_i
! S2(w,x,i,a) += (   -2.00000000) V2(a,a2,x,c0) W7(w,c0,i,a2) 
do s_a2 = 0, nir-1
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a,s_a2) == IEOR(s_x,s_c0) .and. &
IEOR(s_w,s_c0) == IEOR(s_i,s_a2)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- W7(w,c0,i,a2) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i, i_c0, i_a2) =  &
  W7_(s_a2, s_i, s_c0, s_w)%array(i_a2, i_i, i_c0, i_w)
end do
end do
end do
end do
! Z2 <-- V2(a,a2,x,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a2, i_x) =  &
  V2_(s_c0, s_x, s_a2)%array(i_c0, i_x, i_a2)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_w, i_i, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no2_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no3_x0_type1_eri_v &
  (sa, ia, V2, W8, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), W8(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_Xccaa(sleft, W8, nir, nsym, psym) ! -> Xccaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no3_x0_type1_eri_v &
  (sa, ia, h2_i, Xccaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)
deallocate(Xccaa)

end subroutine g_if_sigma_ccov_ccoo_no3_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no3_x0_type1_eri_v &
  (s_a, i_a, V2_, W8_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W8_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_w, i_w, s_c0, i_c0, s_x, i_x, s_i, i_i
! S2(w,x,i,a) += (   -2.00000000) V2(a,a2,w,c0) W8(x,c0,i,a2) 
do s_a2 = 0, nir-1
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_x = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a,s_a2) == IEOR(s_w,s_c0) .and. &
IEOR(s_x,s_c0) == IEOR(s_i,s_a2)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- W8(x,c0,i,a2) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_i, i_c0, i_a2) =  &
  W8_(s_a2, s_i, s_c0, s_x)%array(i_a2, i_i, i_c0, i_x)
end do
end do
end do
end do
! Z2 <-- V2(a,a2,w,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a2, i_w) =  &
  V2_(s_c0, s_w, s_a2)%array(i_c0, i_w, i_a2)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no3_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no4_x0_type1_eri_v &
  (sa, ia, sa1, ia1, V2, W9, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: V2(*), W9(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa1,sa)

call set_symblock_Xaa(sleft, W9, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_ccov_ccoo_no4_x0_type1_eri_v &
  (sa, ia, sa1, ia1, h2_i, Xaa, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no4_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no4_x0_type1_eri_v &
  (s_a, i_a, s_a1, i_a1, V2_, W9_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a1, s_a1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W9_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_a4, i_a4, s_a2, i_a2, s_i, i_i, s_a0, i_a0
! W9(i,a1,a0,a) += (    1.00000000) V2(a,a3,a4,a2) D3(i,a1,a4,a2,a3,a0) 
do s_a3 = 0, nir-1
do s_a4 = 0, nir-1
do s_a2 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a1) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_a3) == IEOR(s_a4,s_a2) .and. &
IEOR(IEOR(s_i,s_a1),s_a4) == IEOR(IEOR(s_a2,s_a3),s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3) > 0) then

! Z1 <-- D3(i,a1,a4,a2,a3,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a4):psym(I_END,I_O, s_a4), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a4 = psym(I_BEGIN, I_O, s_a4), psym(I_END, I_O, s_a4)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a4, i_a2, i_a3) =  &
  D3_(s_a0, s_a3, s_a2, s_a4, s_a1, s_i)%array(i_a0, i_a3, i_a2, i_a4, i_a1, i_i)
end do
end do
end do
end do
end do
! Z2 <-- V2(a,a3,a4,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a4):psym(I_END,I_O, s_a4), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a4 = psym(I_BEGIN, I_O, s_a4), psym(I_END, I_O, s_a4)
Z2_(i_a4, i_a2, i_a3) =  &
  V2_(s_a2, s_a4, s_a3)%array(i_a2, i_a4, i_a3)
end do
end do
end do

! Z3 <-- W9(i,a1,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W9(i,a1,a0,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W9_(s_a0, s_i)%array(i_a0, i_i) = &
    W9_(s_a0, s_i)%array(i_a0, i_i) &
  + Z3_(i_i, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no4_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no4_x1_type1_eri_v &
  (sa, ia, sa1, ia1, T2, W9, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), W9(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa1,sa)

call set_symblock_Xaa(sleft, W9, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no4_x1_type1_eri_v &
  (sa, ia, sa1, ia1, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no4_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no4_x1_type1_eri_v &
  (s_a, i_a, s_a1, i_a1, T2_, W9_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W9_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a0, i_a0, s_i, i_i
! S2(w,x,i,a) += (    2.00000000) T2(x,w,a0,a1) W9(i,a1,a0,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_a0,s_a1) .and. &
IEOR(s_i,s_a1) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W9(i,a1,a0,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0) =  &
  W9_(s_a0, s_i)%array(i_a0, i_i)
end do
end do
! Z2 <-- T2(x,w,a0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x, i_w) =  &
  T2_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x)
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
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
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
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no4_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no5_x0_type1_eri_v &
  (sa, ia, sa0, ia0, V2, W10, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: V2(*), W10(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xaa(sleft, W10, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_ccov_ccoo_no5_x0_type1_eri_v &
  (sa, ia, sa0, ia0, h2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no5_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no5_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, V2_, W10_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W10_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_i, i_i, s_a2, i_a2, s_a1, i_a1
! W10(a1,a0,i,a) += (    1.00000000) V2(a,a3,i,a2) D2(a3,a1,a2,a0) 
do s_a3 = 0, nir-1
do s_i = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(s_a1,s_a0) == IEOR(s_i,s_a) .and. & 
IEOR(s_a,s_a3) == IEOR(s_i,s_a2) .and. &
IEOR(s_a3,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(a,a3,i,a2) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a3, i_a2) =  &
  V2_(s_a2, s_i, s_a3)%array(i_a2, i_i, i_a3)
end do
end do
end do
! Z2 <-- D2(a3,a1,a2,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a2, i_a1) =  &
  D2_(s_a0, s_a2, s_a1, s_a3)%array(i_a0, i_a2, i_a1, i_a3)
end do
end do
end do

! Z3 <-- W10(a1,a0,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! W10(a1,a0,i,a)  <-- Z3
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W10_(s_i, s_a1)%array(i_i, i_a1) = &
    W10_(s_i, s_a1)%array(i_i, i_a1) &
  + Z3_(i_i, i_a1)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no5_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no5_x1_type1_eri_v &
  (sa, ia, sa0, ia0, T2, W10, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), W10(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xaa(sleft, W10, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no5_x1_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_ccov_ccoo_no5_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no5_x1_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, W10_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W10_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a1, i_a1, s_i, i_i
! S2(w,x,i,a) += (    2.00000000) T2(x,w,a1,a0) W10(a1,a0,i,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_a1,s_a0) .and. &
IEOR(s_a1,s_a0) == IEOR(s_i,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- W10(a1,a0,i,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a1) =  &
  W10_(s_i, s_a1)%array(i_i, i_a1)
end do
end do
! Z2 <-- T2(x,w,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_x, i_w) =  &
  T2_(s_a1, s_w, s_x)%array(i_a1, i_w, i_x)
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
                     psym(I_LENGTH,I_O, s_a1),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
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
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no5_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no6_x0_type1_eri_v &
  (sa, ia, sa0, ia0, V2, W11, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: V2(*), W11(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W11, nir, nsym, psym) ! -> Xcc (allocate) 
call g_sigma_ccov_ccoo_no6_x0_type1_eri_v &
  (sa, ia, sa0, ia0, h2_i, Xcc, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no6_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no6_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, V2_, W11_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W11_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_c0, i_c0, s_a1, i_a1
! W11(x,c0,a0,a) += (    1.00000000) V2(a,x,c0,a1) D1(a1,a0) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(s_x,s_c0) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_x) == IEOR(s_c0,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(a,x,c0,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_c0, i_a1) =  &
  V2_(s_a1, s_c0, s_x)%array(i_a1, i_c0, i_x)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do

! Z3 <-- W11(x,c0,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0))

! W11(x,c0,a0,a)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
W11_(s_c0, s_x)%array(i_c0, i_x) = &
    W11_(s_c0, s_x)%array(i_c0, i_x) &
  + Z3_(i_x, i_c0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no6_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no6_x1_type1_eri_v &
  (sa, ia, sa0, ia0, T2, W11, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), W11(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W11, nir, nsym, psym) ! -> Xcc (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no6_x1_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, Xcc, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no6_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no6_x1_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, W11_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W11_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_i, i_i, s_x, i_x
! S2(w,x,i,a) += (   -8.00000000) T2(w,c0,i,a0) W11(x,c0,a0,a) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_w,s_c0) == IEOR(s_i,s_a0) .and. &
IEOR(s_x,s_c0) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(w,c0,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i, i_c0) =  &
  T2_(s_i, s_c0, s_w)%array(i_i, i_c0, i_w)
end do
end do
end do
! Z2 <-- W11(x,c0,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_x) =  &
  W11_(s_c0, s_x)%array(i_c0, i_x)
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_c0),&
                     - 8.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_w, i_i, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no6_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no7_x0_type1_eri_v &
  (sa, ia, sa0, ia0, V2, W12, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: V2(*), W12(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W12, nir, nsym, psym) ! -> Xcc (allocate) 
call g_sigma_ccov_ccoo_no7_x0_type1_eri_v &
  (sa, ia, sa0, ia0, h2_i, Xcc, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no7_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no7_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, V2_, W12_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W12_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_a1, i_a1
! W12(w,c0,a0,a) += (    1.00000000) V2(a,w,c0,a1) D1(a1,a0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(s_w,s_c0) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_c0,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(a,w,c0,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_c0, i_a1) =  &
  V2_(s_a1, s_c0, s_w)%array(i_a1, i_c0, i_w)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do

! Z3 <-- W12(w,c0,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0))

! W12(w,c0,a0,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
W12_(s_c0, s_w)%array(i_c0, i_w) = &
    W12_(s_c0, s_w)%array(i_c0, i_w) &
  + Z3_(i_w, i_c0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no7_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no7_x1_type1_eri_v &
  (sa, ia, sa0, ia0, T2, W12, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), W12(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W12, nir, nsym, psym) ! -> Xcc (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no7_x1_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, Xcc, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no7_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no7_x1_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, W12_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W12_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_c0, i_c0, s_i, i_i, s_w, i_w
! S2(w,x,i,a) += (    4.00000000) T2(x,c0,i,a0) W12(w,c0,a0,a) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_c0) == IEOR(s_i,s_a0) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(x,c0,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_i, i_c0) =  &
  T2_(s_i, s_c0, s_x)%array(i_i, i_c0, i_x)
end do
end do
end do
! Z2 <-- W12(w,c0,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  W12_(s_c0, s_w)%array(i_c0, i_w)
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_ccov_ccoo_no7_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no8_x0_type1_eri_v &
  (sa, ia, V2, W13, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), W13(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_Xcca(sleft, W13, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_ccov_ccoo_no8_x0_type1_eri_v &
  (sa, ia, h2_i, Xcca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no8_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no8_x0_type1_eri_v &
  (s_a, i_a, V2_, W13_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W13_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_a1, i_a1, s_a0, i_a0
! W13(w,c0,a0,a) += (    1.00000000) V2(a,w,c0,a1) D1(a1,a0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_c0) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_c0,s_a1) .and. &
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
! Z2 <-- V2(a,w,c0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w, i_c0) =  &
  V2_(s_a1, s_c0, s_w)%array(i_a1, i_c0, i_w)
end do
end do
end do

! Z3 <-- W13(w,c0,a0,a) 
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

! W13(w,c0,a0,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W13_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w) = &
    W13_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w) &
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

end subroutine g_sigma_ccov_ccoo_no8_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no8_x1_type1_eri_v &
  (sa, ia, si, ii, T2, W13, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, si, ii
real(kind=8), intent(inout) :: T2(*), W13(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(si, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcca(sleft, W13, nir, nsym, psym) ! -> Xcca (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no8_x1_type1_eri_v &
  (sa, ia, si, ii, av2_i, Xcca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no8_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no8_x1_type1_eri_v &
  (s_a, i_a, s_i, i_i, T2_, W13_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_a, s_a
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

integer :: s_x, i_x, s_c0, i_c0, s_a0, i_a0, s_w, i_w
! S2(w,x,i,a) += (   -2.00000000) T2(x,c0,a0,i) W13(w,c0,a0,a) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_c0) == IEOR(s_a0,s_i) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(x,c0,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_c0, i_a0) =  &
  T2_(s_a0, s_c0, s_x)%array(i_a0, i_c0, i_x)
end do
end do
end do
! Z2 <-- W13(w,c0,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0, i_w) =  &
  W13_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x) * &
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

end subroutine g_sigma_ccov_ccoo_no8_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no9_x0_type1_eri_v &
  (sa, ia, V2, W14, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), W14(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_Xcca(sleft, W14, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_ccov_ccoo_no9_x0_type1_eri_v &
  (sa, ia, h2_i, Xcca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no9_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no9_x0_type1_eri_v &
  (s_a, i_a, V2_, W14_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W14_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_c0, i_c0, s_a1, i_a1, s_a0, i_a0
! W14(x,c0,a0,a) += (    1.00000000) V2(a,x,c0,a1) D1(a1,a0) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_c0) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_x) == IEOR(s_c0,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0) > 0 .and. &
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
! Z2 <-- V2(a,x,c0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_x, i_c0) =  &
  V2_(s_a1, s_c0, s_x)%array(i_a1, i_c0, i_x)
end do
end do
end do

! Z3 <-- W14(x,c0,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W14(x,c0,a0,a)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W14_(s_a0, s_c0, s_x)%array(i_a0, i_c0, i_x) = &
    W14_(s_a0, s_c0, s_x)%array(i_a0, i_c0, i_x) &
  + Z3_(i_a0, i_x, i_c0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no9_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no9_x1_type1_eri_v &
  (sa, ia, si, ii, T2, W14, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, si, ii
real(kind=8), intent(inout) :: T2(*), W14(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(si, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcca(sleft, W14, nir, nsym, psym) ! -> Xcca (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no9_x1_type1_eri_v &
  (sa, ia, si, ii, av2_i, Xcca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no9_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no9_x1_type1_eri_v &
  (s_a, i_a, s_i, i_i, T2_, W14_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_a, s_a
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

integer :: s_w, i_w, s_c0, i_c0, s_a0, i_a0, s_x, i_x
! S2(w,x,i,a) += (    4.00000000) T2(w,c0,a0,i) W14(x,c0,a0,a) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_w,s_c0) == IEOR(s_a0,s_i) .and. &
IEOR(s_x,s_c0) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W14(x,c0,a0,a) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_c0, i_a0) =  &
  W14_(s_a0, s_c0, s_x)%array(i_a0, i_c0, i_x)
end do
end do
end do
! Z2 <-- T2(w,c0,a0,i) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0, i_w) =  &
  T2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x) * &
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

end subroutine g_sigma_ccov_ccoo_no9_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no10_x0_type1_eri_v &
  (sa, ia, sa0, ia0, V2, W15, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: V2(*), W15(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W15, nir, nsym, psym) ! -> Xcc (allocate) 
call g_sigma_ccov_ccoo_no10_x0_type1_eri_v &
  (sa, ia, sa0, ia0, h2_i, Xcc, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no10_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no10_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, V2_, W15_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W15_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_w, i_w, s_c0, i_c0
! W15(w,c0,a0,a) += (    1.00000000) V2(a,a1,w,c0) D1(a1,a0) 
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_c0 = 0, nir-1
if( &
IEOR(s_w,s_c0) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_a1) == IEOR(s_w,s_c0) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(a,a1,w,c0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_c0, i_a1) =  &
  V2_(s_c0, s_w, s_a1)%array(i_c0, i_w, i_a1)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do

! Z3 <-- W15(w,c0,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0))

! W15(w,c0,a0,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
W15_(s_c0, s_w)%array(i_c0, i_w) = &
    W15_(s_c0, s_w)%array(i_c0, i_w) &
  + Z3_(i_w, i_c0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_c0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no10_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no10_x1_type1_eri_v &
  (sa, ia, sa0, ia0, T2, W15, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), W15(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W15, nir, nsym, psym) ! -> Xcc (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no10_x1_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, Xcc, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no10_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no10_x1_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, W15_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W15_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_c0, i_c0, s_i, i_i, s_w, i_w
! S2(w,x,i,a) += (   -2.00000000) T2(x,c0,i,a0) W15(w,c0,a0,a) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_c0) == IEOR(s_i,s_a0) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(x,c0,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_i, i_c0) =  &
  T2_(s_i, s_c0, s_x)%array(i_i, i_c0, i_x)
end do
end do
end do
! Z2 <-- W15(w,c0,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  W15_(s_c0, s_w)%array(i_c0, i_w)
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_ccov_ccoo_no10_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no11_x0_type1_eri_v &
  (sa, ia, sa0, ia0, V2, W16, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: V2(*), W16(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W16, nir, nsym, psym) ! -> Xcc (allocate) 
call g_sigma_ccov_ccoo_no11_x0_type1_eri_v &
  (sa, ia, sa0, ia0, h2_i, Xcc, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no11_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no11_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, V2_, W16_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W16_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_x, i_x, s_c0, i_c0
! W16(x,c0,a0,a) += (    1.00000000) V2(a,a1,x,c0) D1(a1,a0) 
do s_a1 = 0, nir-1
do s_x = 0, nir-1
do s_c0 = 0, nir-1
if( &
IEOR(s_x,s_c0) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_a1) == IEOR(s_x,s_c0) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(a,a1,x,c0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_c0, i_a1) =  &
  V2_(s_c0, s_x, s_a1)%array(i_c0, i_x, i_a1)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do

! Z3 <-- W16(x,c0,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0))

! W16(x,c0,a0,a)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
W16_(s_c0, s_x)%array(i_c0, i_x) = &
    W16_(s_c0, s_x)%array(i_c0, i_x) &
  + Z3_(i_x, i_c0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no11_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no11_x1_type1_eri_v &
  (sa, ia, sa0, ia0, T2, W16, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), W16(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W16, nir, nsym, psym) ! -> Xcc (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no11_x1_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, Xcc, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no11_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no11_x1_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, W16_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
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

integer :: s_w, i_w, s_c0, i_c0, s_i, i_i, s_x, i_x
! S2(w,x,i,a) += (    4.00000000) T2(w,c0,i,a0) W16(x,c0,a0,a) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_w,s_c0) == IEOR(s_i,s_a0) .and. &
IEOR(s_x,s_c0) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(w,c0,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i, i_c0) =  &
  T2_(s_i, s_c0, s_w)%array(i_i, i_c0, i_w)
end do
end do
end do
! Z2 <-- W16(x,c0,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_x) =  &
  W16_(s_c0, s_x)%array(i_c0, i_x)
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_c0),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_w, i_i, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no11_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no12_x0_type1_eri_v &
  (sa, ia, V2, W17, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), W17(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_Xcca(sleft, W17, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_ccov_ccoo_no12_x0_type1_eri_v &
  (sa, ia, h2_i, Xcca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no12_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no12_x0_type1_eri_v &
  (s_a, i_a, V2_, W17_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W17_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_w, i_w, s_c0, i_c0, s_a0, i_a0
! W17(w,c0,a0,a) += (    1.00000000) V2(a,a1,w,c0) D1(a1,a0) 
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_c0) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_a1) == IEOR(s_w,s_c0) .and. &
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
! Z2 <-- V2(a,a1,w,c0) 
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

! Z3 <-- W17(w,c0,a0,a) 
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

! W17(w,c0,a0,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W17_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w) = &
    W17_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w) &
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

end subroutine g_sigma_ccov_ccoo_no12_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no12_x1_type1_eri_v &
  (sa, ia, si, ii, T2, W17, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, si, ii
real(kind=8), intent(inout) :: T2(*), W17(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(si, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcca(sleft, W17, nir, nsym, psym) ! -> Xcca (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no12_x1_type1_eri_v &
  (sa, ia, si, ii, av2_i, Xcca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no12_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no12_x1_type1_eri_v &
  (s_a, i_a, s_i, i_i, T2_, W17_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W17_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_c0, i_c0, s_a0, i_a0, s_w, i_w
! S2(w,x,i,a) += (    4.00000000) T2(x,c0,a0,i) W17(w,c0,a0,a) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_c0) == IEOR(s_a0,s_i) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(x,c0,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_c0, i_a0) =  &
  T2_(s_a0, s_c0, s_x)%array(i_a0, i_c0, i_x)
end do
end do
end do
! Z2 <-- W17(w,c0,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0, i_w) =  &
  W17_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x) * &
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

end subroutine g_sigma_ccov_ccoo_no12_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no13_x0_type1_eri_v &
  (sa, ia, V2, W18, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), W18(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_Xcca(sleft, W18, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_ccov_ccoo_no13_x0_type1_eri_v &
  (sa, ia, h2_i, Xcca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no13_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no13_x0_type1_eri_v &
  (s_a, i_a, V2_, W18_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W18_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_x, i_x, s_c0, i_c0, s_a0, i_a0
! W18(x,c0,a0,a) += (    1.00000000) V2(a,a1,x,c0) D1(a1,a0) 
do s_a1 = 0, nir-1
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_c0) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_a1) == IEOR(s_x,s_c0) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0) > 0 .and. &
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
! Z2 <-- V2(a,a1,x,c0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_x, i_c0) =  &
  V2_(s_c0, s_x, s_a1)%array(i_c0, i_x, i_a1)
end do
end do
end do

! Z3 <-- W18(x,c0,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W18(x,c0,a0,a)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W18_(s_a0, s_c0, s_x)%array(i_a0, i_c0, i_x) = &
    W18_(s_a0, s_c0, s_x)%array(i_a0, i_c0, i_x) &
  + Z3_(i_a0, i_x, i_c0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no13_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no13_x1_type1_eri_v &
  (sa, ia, si, ii, T2, W18, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, si, ii
real(kind=8), intent(inout) :: T2(*), W18(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(si, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcca(sleft, W18, nir, nsym, psym) ! -> Xcca (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no13_x1_type1_eri_v &
  (sa, ia, si, ii, av2_i, Xcca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no13_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no13_x1_type1_eri_v &
  (s_a, i_a, s_i, i_i, T2_, W18_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_a, s_a
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

integer :: s_w, i_w, s_c0, i_c0, s_a0, i_a0, s_x, i_x
! S2(w,x,i,a) += (   -2.00000000) T2(w,c0,a0,i) W18(x,c0,a0,a) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_w,s_c0) == IEOR(s_a0,s_i) .and. &
IEOR(s_x,s_c0) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W18(x,c0,a0,a) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_c0, i_a0) =  &
  W18_(s_a0, s_c0, s_x)%array(i_a0, i_c0, i_x)
end do
end do
end do
! Z2 <-- T2(w,c0,a0,i) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0, i_w) =  &
  T2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x) * &
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

end subroutine g_sigma_ccov_ccoo_no13_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no14_x0_type1_eri_v &
  (sa, ia, sa0, ia0, V2, W19, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: V2(*), W19
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa0,sa)

call g_sigma_ccov_ccoo_no14_x0_type1_eri_v &
  (sa, ia, sa0, ia0, h2_i, W19, d2, nir, nsym, psym, flops)

deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no14_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no14_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, V2_, W19_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W19_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a3, i_a3, s_a1, i_a1
! W19(a0,a) += (    1.00000000) V2(a,a2,a3,a1) D2(a3,a1,a2,a0) 
do s_a2 = 0, nir-1
do s_a3 = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(s_a0,s_a) == 0 .and. & 
IEOR(s_a,s_a2) == IEOR(s_a3,s_a1) .and. &
IEOR(s_a3,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(a,a2,a3,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_a3, i_a1) =  &
  V2_(s_a1, s_a3, s_a2)%array(i_a1, i_a3, i_a2)
end do
end do
end do
! Z2 <-- D2(a3,a1,a2,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a3, i_a1) =  &
  D2_(s_a0, s_a2, s_a1, s_a3)%array(i_a0, i_a2, i_a1, i_a3)
end do
end do
end do

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', 1,&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     1,&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     1)

! W19(a0,a)  <-- Z3
W19_ = &
    W19_ &
  + Z3_

! Flop count
flops = flops + 1 * &
                1 * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no14_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no14_x1_type1_eri_v &
  (sa, ia, sa0, ia0, T2, W19, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), W19, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no14_x1_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, W19, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_ccov_ccoo_no14_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no14_x1_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, W19_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W19_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_i, i_i
! S2(w,x,i,a) += (    2.00000000) T2(x,w,i,a0) W19(a0,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_i,s_a0) .and. &
IEOR(s_a0,s_a) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0) then

! Z1 <-- T2(x,w,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_w, i_i) =  &
  T2_(s_i, s_w, s_x)%array(i_i, i_w, i_x)
end do
end do
end do
! Z2 <-- W19(a0,a) 
Z2_ =  &
  W19_

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     1,&
                     2.00000000d+00, &
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

end subroutine g_sigma_ccov_ccoo_no14_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no15_x0_type1_eri_v &
  (sa, ia, V2, W20, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), W20(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_Xa(sleft, W20, nir, nsym, psym) ! -> Xa (allocate) 
call g_sigma_ccov_ccoo_no15_x0_type1_eri_v &
  (sa, ia, h2_i, Xa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xa)

end subroutine g_if_sigma_ccov_ccoo_no15_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no15_x0_type1_eri_v &
  (s_a, i_a, V2_, W20_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W20_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a3, i_a3, s_a1, i_a1, s_a0, i_a0
! W20(a0,a) += (    1.00000000) V2(a,a2,a3,a1) D2(a3,a1,a2,a0) 
do s_a2 = 0, nir-1
do s_a3 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_a0,s_a) == 0 .and. & 
IEOR(s_a,s_a2) == IEOR(s_a3,s_a1) .and. &
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
! Z2 <-- V2(a,a2,a3,a1) 
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

! Z3 <-- W20(a0,a) 
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

! W20(a0,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W20_(s_a0)%array(i_a0) = &
    W20_(s_a0)%array(i_a0) &
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

end subroutine g_sigma_ccov_ccoo_no15_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no15_x1_type1_eri_v &
  (sa, ia, si, ii, T2, W20, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, si, ii
real(kind=8), intent(inout) :: T2(*), W20(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(si, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xa(sleft, W20, nir, nsym, psym) ! -> Xa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no15_x1_type1_eri_v &
  (sa, ia, si, ii, av2_i, Xa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xa)

end subroutine g_if_sigma_ccov_ccoo_no15_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no15_x1_type1_eri_v &
  (s_a, i_a, s_i, i_i, T2_, W20_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W20_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_a0, i_a0
! S2(w,x,i,a) += (   -4.00000000) T2(x,w,a0,i) W20(a0,a) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_x,s_w) == IEOR(s_a0,s_i) .and. &
IEOR(s_a0,s_a) == 0) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(x,w,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_w, i_a0) =  &
  T2_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x)
end do
end do
end do
! Z2 <-- W20(a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0) =  &
  W20_(s_a0)%array(i_a0)
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w),&
                     1,&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) * &
                1 * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no15_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no16_x0_type1_eri_v &
  (sa, ia, sa0, ia0, T2, V2, W21, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), W21(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W21, nir, nsym, psym) ! -> Xcc (allocate) 
call g_sigma_ccov_ccoo_no16_x0_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, h2_i, Xcc, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no16_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no16_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, W21_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W21_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_c0, i_c0, s_a1, i_a1, s_w, i_w
! W21(w,x,a0,a) += (    1.00000000) V2(a,x,c0,a1) T2(w,c0,a1,a0) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_x) == IEOR(s_c0,s_a1) .and. &
IEOR(s_w,s_c0) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(a,x,c0,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_c0, i_a1) =  &
  V2_(s_a1, s_c0, s_x)%array(i_a1, i_c0, i_x)
end do
end do
end do
! Z2 <-- T2(w,c0,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a1, i_w) =  &
  T2_(s_a1, s_c0, s_w)%array(i_a1, i_c0, i_w)
end do
end do
end do

! Z3 <-- W21(w,x,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x))

! W21(w,x,a0,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
W21_(s_x, s_w)%array(i_x, i_w) = &
    W21_(s_x, s_w)%array(i_x, i_w) &
  + Z3_(i_x, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no16_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no16_x1_type1_eri_v &
  (sa, ia, sa0, ia0, W21, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: W21(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W21, nir, nsym, psym) ! -> Xcc (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no16_x1_type1_eri_v &
  (sa, ia, sa0, ia0, Xcc, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no16_x1_type1_eri_v



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
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ccoo_no16_x1_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, W21_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W21_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_w, i_w, s_x, i_x
! S2(w,x,i,a) += (    4.00000000) D1(i,a0) W21(w,x,a0,a) 
do s_i = 0, nir-1
do s_w = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_x) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0) then

! Z1 <-- D1(i,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
! Z2 <-- W21(w,x,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z2_(i_w, i_x) =  &
  W21_(s_x, s_w)%array(i_x, i_w)
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     1,&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_i, i_w, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no16_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no17_x0_type1_eri_v &
  (sa, ia, sa0, ia0, T2, V2, W22, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), W22(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W22, nir, nsym, psym) ! -> Xcc (allocate) 
call g_sigma_ccov_ccoo_no17_x0_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, h2_i, Xcc, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no17_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no17_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, W22_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a0, s_a0
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

integer :: s_w, i_w, s_c0, i_c0, s_a1, i_a1, s_x, i_x
! W22(x,w,a0,a) += (    1.00000000) V2(a,w,c0,a1) T2(x,c0,a1,a0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_x,s_w) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_c0,s_a1) .and. &
IEOR(s_x,s_c0) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(a,w,c0,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_c0, i_a1) =  &
  V2_(s_a1, s_c0, s_w)%array(i_a1, i_c0, i_w)
end do
end do
end do
! Z2 <-- T2(x,c0,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a1, i_x) =  &
  T2_(s_a1, s_c0, s_x)%array(i_a1, i_c0, i_x)
end do
end do
end do

! Z3 <-- W22(x,w,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w))

! W22(x,w,a0,a)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
W22_(s_w, s_x)%array(i_w, i_x) = &
    W22_(s_w, s_x)%array(i_w, i_x) &
  + Z3_(i_w, i_x)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no17_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no17_x1_type1_eri_v &
  (sa, ia, sa0, ia0, W22, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: W22(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W22, nir, nsym, psym) ! -> Xcc (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no17_x1_type1_eri_v &
  (sa, ia, sa0, ia0, Xcc, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no17_x1_type1_eri_v



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
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ccoo_no17_x1_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, W22_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W22_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_x, i_x, s_w, i_w
! S2(w,x,i,a) += (   -2.00000000) D1(i,a0) W22(x,w,a0,a) 
do s_i = 0, nir-1
do s_x = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_x,s_w) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0) then

! Z1 <-- D1(i,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
! Z2 <-- W22(x,w,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z2_(i_x, i_w) =  &
  W22_(s_w, s_x)%array(i_w, i_x)
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w),&
                     1,&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
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
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no17_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no18_x0_type1_eri_v &
  (sa, ia, sa1, ia1, T2, V2, W23, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), V2(*), W23(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcca(sleft, W23, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_ccov_ccoo_no18_x0_type1_eri_v &
  (sa, ia, sa1, ia1, av2_i, h2_i, Xcca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no18_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no18_x0_type1_eri_v &
  (s_a, i_a, s_a1, i_a1, T2_, V2_, W23_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a1, s_a1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W23_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_x, i_x, s_a0, i_a0
! W23(x,w,a0,a) += (    1.00000000) V2(a,w,c0,a1) T2(x,c0,a0,a1) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_x = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_w) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_c0,s_a1) .and. &
IEOR(s_x,s_c0) == IEOR(s_a0,s_a1)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(x,c0,a0,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_a0, i_c0) =  &
  T2_(s_a0, s_c0, s_x)%array(i_a0, i_c0, i_x)
end do
end do
end do
! Z2 <-- V2(a,w,c0,a1) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  V2_(s_a1, s_c0, s_w)%array(i_a1, i_c0, i_w)
end do
end do

! Z3 <-- W23(x,w,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_a0))

! W23(x,w,a0,a)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W23_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x) = &
    W23_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x) &
  + Z3_(i_x, i_a0, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_a0) * &
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

end subroutine g_sigma_ccov_ccoo_no18_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no18_x1_type1_eri_v &
  (sa, ia, W23, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: W23(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa

call set_symblock_Xcca(sleft, W23, nir, nsym, psym) ! -> Xcca (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no18_x1_type1_eri_v &
  (sa, ia, Xcca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no18_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no18_x1_type1_eri_v &
  (s_a, i_a, W23_, S2_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W23_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_x, i_x, s_w, i_w
! S2(w,x,i,a) += (    4.00000000) D1(i,a0) W23(x,w,a0,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_x = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_x,s_w) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
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
! Z2 <-- W23(x,w,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x, i_w) =  &
  W23_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x)
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
                     psym(I_LENGTH,I_O, s_a0),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
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
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no18_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no19_x0_type1_eri_v &
  (sa, ia, sa1, ia1, T2, V2, W24, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), V2(*), W24(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcca(sleft, W24, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_ccov_ccoo_no19_x0_type1_eri_v &
  (sa, ia, sa1, ia1, av2_i, h2_i, Xcca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no19_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no19_x0_type1_eri_v &
  (s_a, i_a, s_a1, i_a1, T2_, V2_, W24_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a1, s_a1
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

integer :: s_x, i_x, s_c0, i_c0, s_w, i_w, s_a0, i_a0
! W24(w,x,a0,a) += (    1.00000000) V2(a,x,c0,a1) T2(w,c0,a0,a1) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_x) == IEOR(s_c0,s_a1) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_a1)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(w,c0,a0,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a0, i_c0) =  &
  T2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do
! Z2 <-- V2(a,x,c0,a1) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_x) =  &
  V2_(s_a1, s_c0, s_x)%array(i_a1, i_c0, i_x)
end do
end do

! Z3 <-- W24(w,x,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0))

! W24(w,x,a0,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W24_(s_a0, s_x, s_w)%array(i_a0, i_x, i_w) = &
    W24_(s_a0, s_x, s_w)%array(i_a0, i_x, i_w) &
  + Z3_(i_w, i_a0, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no19_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no19_x1_type1_eri_v &
  (sa, ia, W24, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: W24(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa

call set_symblock_Xcca(sleft, W24, nir, nsym, psym) ! -> Xcca (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no19_x1_type1_eri_v &
  (sa, ia, Xcca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no19_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no19_x1_type1_eri_v &
  (s_a, i_a, W24_, S2_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W24_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_w, i_w, s_x, i_x
! S2(w,x,i,a) += (   -8.00000000) D1(i,a0) W24(w,x,a0,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_x) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0 .and. &
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
! Z2 <-- W24(w,x,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w, i_x) =  &
  W24_(s_a0, s_x, s_w)%array(i_a0, i_x, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 8.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_i, i_w, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no19_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no20_x0_type1_eri_v &
  (sa, ia, sa0, ia0, T2, V2, W27, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), W27(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W27, nir, nsym, psym) ! -> Xcc (allocate) 
call g_sigma_ccov_ccoo_no20_x0_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, h2_i, Xcc, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no20_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no20_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, W27_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W27_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_w, i_w, s_c0, i_c0, s_x, i_x
! W27(x,w,a0,a) += (    1.00000000) V2(a,a1,w,c0) T2(x,c0,a1,a0) 
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_x,s_w) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_a1) == IEOR(s_w,s_c0) .and. &
IEOR(s_x,s_c0) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(a,a1,w,c0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a1, i_c0) =  &
  V2_(s_c0, s_w, s_a1)%array(i_c0, i_w, i_a1)
end do
end do
end do
! Z2 <-- T2(x,c0,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_c0, i_x) =  &
  T2_(s_a1, s_c0, s_x)%array(i_a1, i_c0, i_x)
end do
end do
end do

! Z3 <-- W27(x,w,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w))

! W27(x,w,a0,a)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
W27_(s_w, s_x)%array(i_w, i_x) = &
    W27_(s_w, s_x)%array(i_w, i_x) &
  + Z3_(i_w, i_x)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no20_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no20_x1_type1_eri_v &
  (sa, ia, sa0, ia0, W27, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: W27(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W27, nir, nsym, psym) ! -> Xcc (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no20_x1_type1_eri_v &
  (sa, ia, sa0, ia0, Xcc, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no20_x1_type1_eri_v



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
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ccoo_no20_x1_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, W27_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W27_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_x, i_x, s_w, i_w
! S2(w,x,i,a) += (    4.00000000) D1(i,a0) W27(x,w,a0,a) 
do s_i = 0, nir-1
do s_x = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_x,s_w) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0) then

! Z1 <-- D1(i,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
! Z2 <-- W27(x,w,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z2_(i_x, i_w) =  &
  W27_(s_w, s_x)%array(i_w, i_x)
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w),&
                     1,&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
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
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no20_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no21_x0_type1_eri_v &
  (sa, ia, sa0, ia0, T2, V2, W28, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), W28(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W28, nir, nsym, psym) ! -> Xcc (allocate) 
call g_sigma_ccov_ccoo_no21_x0_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, h2_i, Xcc, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no21_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no21_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, W28_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a0, s_a0
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

integer :: s_a1, i_a1, s_x, i_x, s_c0, i_c0, s_w, i_w
! W28(w,x,a0,a) += (    1.00000000) V2(a,a1,x,c0) T2(w,c0,a1,a0) 
do s_a1 = 0, nir-1
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_a1) == IEOR(s_x,s_c0) .and. &
IEOR(s_w,s_c0) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(a,a1,x,c0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_a1, i_c0) =  &
  V2_(s_c0, s_x, s_a1)%array(i_c0, i_x, i_a1)
end do
end do
end do
! Z2 <-- T2(w,c0,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_c0, i_w) =  &
  T2_(s_a1, s_c0, s_w)%array(i_a1, i_c0, i_w)
end do
end do
end do

! Z3 <-- W28(w,x,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x))

! W28(w,x,a0,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
W28_(s_x, s_w)%array(i_x, i_w) = &
    W28_(s_x, s_w)%array(i_x, i_w) &
  + Z3_(i_x, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no21_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no21_x1_type1_eri_v &
  (sa, ia, sa0, ia0, W28, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: W28(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa0,sa)

call set_symblock_Xcc(sleft, W28, nir, nsym, psym) ! -> Xcc (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no21_x1_type1_eri_v &
  (sa, ia, sa0, ia0, Xcc, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcc)

end subroutine g_if_sigma_ccov_ccoo_no21_x1_type1_eri_v



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
!    -- allRDM.first false
!    -- precedence   -1
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ccov_ccoo_no21_x1_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, W28_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W28_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_w, i_w, s_x, i_x
! S2(w,x,i,a) += (   -2.00000000) D1(i,a0) W28(w,x,a0,a) 
do s_i = 0, nir-1
do s_w = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_x) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0) then

! Z1 <-- D1(i,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i) =  &
  D1_(s_a0, s_i)%array(i_a0, i_i)
end do
! Z2 <-- W28(w,x,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z2_(i_w, i_x) =  &
  W28_(s_x, s_w)%array(i_x, i_w)
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     1,&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_i, i_w, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no21_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no22_x0_type1_eri_v &
  (sa, ia, sa1, ia1, T2, V2, W29, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), V2(*), W29(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcca(sleft, W29, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_ccov_ccoo_no22_x0_type1_eri_v &
  (sa, ia, sa1, ia1, av2_i, h2_i, Xcca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no22_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no22_x0_type1_eri_v &
  (s_a, i_a, s_a1, i_a1, T2_, V2_, W29_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a1, s_a1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W29_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_x, i_x, s_a0, i_a0
! W29(x,w,a0,a) += (    1.00000000) V2(a,a1,w,c0) T2(x,c0,a0,a1) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_x = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_x,s_w) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_a1) == IEOR(s_w,s_c0) .and. &
IEOR(s_x,s_c0) == IEOR(s_a0,s_a1)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(x,c0,a0,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_a0, i_c0) =  &
  T2_(s_a0, s_c0, s_x)%array(i_a0, i_c0, i_x)
end do
end do
end do
! Z2 <-- V2(a,a1,w,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  V2_(s_c0, s_w, s_a1)%array(i_c0, i_w, i_a1)
end do
end do

! Z3 <-- W29(x,w,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_a0))

! W29(x,w,a0,a)  <-- Z3
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W29_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x) = &
    W29_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x) &
  + Z3_(i_x, i_a0, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_a0) * &
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

end subroutine g_sigma_ccov_ccoo_no22_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no22_x1_type1_eri_v &
  (sa, ia, W29, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: W29(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa

call set_symblock_Xcca(sleft, W29, nir, nsym, psym) ! -> Xcca (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no22_x1_type1_eri_v &
  (sa, ia, Xcca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no22_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no22_x1_type1_eri_v &
  (s_a, i_a, W29_, S2_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W29_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_x, i_x, s_w, i_w
! S2(w,x,i,a) += (   -2.00000000) D1(i,a0) W29(x,w,a0,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_x = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_x,s_w) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
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
! Z2 <-- W29(x,w,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_x, i_w) =  &
  W29_(s_a0, s_w, s_x)%array(i_a0, i_w, i_x)
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
                     psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
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
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no22_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no23_x0_type1_eri_v &
  (sa, ia, sa1, ia1, T2, V2, W30, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), V2(*), W30(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcca(sleft, W30, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_ccov_ccoo_no23_x0_type1_eri_v &
  (sa, ia, sa1, ia1, av2_i, h2_i, Xcca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no23_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no23_x0_type1_eri_v &
  (s_a, i_a, s_a1, i_a1, T2_, V2_, W30_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a1, s_a1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W30_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_c0, i_c0, s_w, i_w, s_a0, i_a0
! W30(w,x,a0,a) += (    1.00000000) V2(a,a1,x,c0) T2(w,c0,a0,a1) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_a1) == IEOR(s_x,s_c0) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_a1)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(w,c0,a0,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a0, i_c0) =  &
  T2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do
! Z2 <-- V2(a,a1,x,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_x) =  &
  V2_(s_c0, s_x, s_a1)%array(i_c0, i_x, i_a1)
end do
end do

! Z3 <-- W30(w,x,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0))

! W30(w,x,a0,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W30_(s_a0, s_x, s_w)%array(i_a0, i_x, i_w) = &
    W30_(s_a0, s_x, s_w)%array(i_a0, i_x, i_w) &
  + Z3_(i_w, i_a0, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no23_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no23_x1_type1_eri_v &
  (sa, ia, W30, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: W30(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa

call set_symblock_Xcca(sleft, W30, nir, nsym, psym) ! -> Xcca (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no23_x1_type1_eri_v &
  (sa, ia, Xcca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcca)

end subroutine g_if_sigma_ccov_ccoo_no23_x1_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no23_x1_type1_eri_v &
  (s_a, i_a, W30_, S2_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W30_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_w, i_w, s_x, i_x
! S2(w,x,i,a) += (    4.00000000) D1(i,a0) W30(w,x,a0,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_x) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0 .and. &
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
! Z2 <-- W30(w,x,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w, i_x) =  &
  W30_(s_a0, s_x, s_w)%array(i_a0, i_x, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_O, s_a0),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_i, i_w, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no23_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no24_x0_type1_eri_v &
  (sa, ia, sa0, ia0, T2, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no24_x0_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, h2_i, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no24_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no24_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_c0, i_c0, s_w, i_w, s_i, i_i
! S2(w,x,i,a) += (   16.00000000) V2(a,x,c0,a0) T2(w,c0,i,a0) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a,s_x) == IEOR(s_c0,s_a0) .and. &
IEOR(s_w,s_c0) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(w,c0,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i, i_c0) =  &
  T2_(s_i, s_c0, s_w)%array(i_i, i_c0, i_w)
end do
end do
end do
! Z2 <-- V2(a,x,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_x) =  &
  V2_(s_a0, s_c0, s_x)%array(i_a0, i_c0, i_x)
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_c0),&
                     16.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_w, i_i, i_x)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_x) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccov_ccoo_no24_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no25_x0_type1_eri_v &
  (sa, ia, sa0, ia0, T2, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no25_x0_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, h2_i, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no25_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no25_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_x, i_x, s_i, i_i
! S2(w,x,i,a) += (   -8.00000000) V2(a,w,c0,a0) T2(x,c0,i,a0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_x = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_c0,s_a0) .and. &
IEOR(s_x,s_c0) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(x,c0,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_i, i_c0) =  &
  T2_(s_i, s_c0, s_x)%array(i_i, i_c0, i_x)
end do
end do
end do
! Z2 <-- V2(a,w,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  V2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     - 8.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_ccov_ccoo_no25_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no26_x0_type1_eri_v &
  (sa, ia, si, ii, T2, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, si, ii
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(si, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no26_x0_type1_eri_v &
  (sa, ia, si, ii, av2_i, h2_i, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no26_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no26_x0_type1_eri_v &
  (s_a, i_a, s_i, i_i, T2_, V2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_i, s_i
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_a0, i_a0, s_x, i_x
! S2(w,x,i,a) += (    4.00000000) V2(a,w,c0,a0) T2(x,c0,a0,i) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_c0,s_a0) .and. &
IEOR(s_x,s_c0) == IEOR(s_a0,s_i)) then

if(psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(x,c0,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_c0, i_a0) =  &
  T2_(s_a0, s_c0, s_x)%array(i_a0, i_c0, i_x)
end do
end do
end do
! Z2 <-- V2(a,w,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0, i_w) =  &
  V2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x) * &
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

end subroutine g_sigma_ccov_ccoo_no26_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccov_ccoo_no27_x0_type1_eri_v &
  (sa, ia, si, ii, T2, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, si, ii
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(si, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccov_ccoo_no27_x0_type1_eri_v &
  (sa, ia, si, ii, av2_i, h2_i, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_ccov_ccoo_no27_x0_type1_eri_v



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
subroutine g_sigma_ccov_ccoo_no27_x0_type1_eri_v &
  (s_a, i_a, s_i, i_i, T2_, V2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_i, s_i
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_c0, i_c0, s_a0, i_a0, s_w, i_w
! S2(w,x,i,a) += (   -8.00000000) V2(a,x,c0,a0) T2(w,c0,a0,i) 
do s_x = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_a) .and. & 
IEOR(s_a,s_x) == IEOR(s_c0,s_a0) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_i)) then

if(psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(a,x,c0,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
Z1_(i_x, i_c0, i_a0) =  &
  V2_(s_a0, s_c0, s_x)%array(i_a0, i_c0, i_x)
end do
end do
end do
! Z2 <-- T2(w,c0,a0,i) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0, i_w) =  &
  T2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     - 8.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_x),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_x))

! S2(w,x,i,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) = &
    S2_(s_i, s_x, s_w)%array(i_i, i_x, i_w) &
  + Z3_(i_x, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_x) * &
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

end subroutine g_sigma_ccov_ccoo_no27_x0_type1_eri_v

