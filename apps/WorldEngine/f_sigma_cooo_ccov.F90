#include <sci/icmr/fsrc/f_mr.fh>


!  ___________                __               
!  \_   _____/____    _____ _/  |_  ____      
!   |    __)_/ __ \  /     \\   __\/  _ \ 
!   |     \ \  ___/ |  Y Y  \|  | (  <_> )  
!   \___  /  \___  >|__|_|  /|__|  \____/   
!       \/       \/       \/                

!                                    Generated date : Sun Apr 20 10:26:17 2014



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no0_x0_type0_noeri &
  (sa0, ia0, T2, W0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0
real(kind=8), intent(inout) :: T2(*), W0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa0

call set_symblock_Xc(sleft, W0, nir, nsym, psym) ! -> Xc (allocate) 
call g_sigma_cooo_ccov_no0_x0_type0_noeri &
  (sa0, ia0, av2_i, Xc, fc1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no0_x0_type0_noeri



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
subroutine g_sigma_cooo_ccov_no0_x0_type0_noeri &
  (s_a0, i_a0, T2_, W0_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W0_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_v0, i_v0
! W0(w,a0) += (    1.00000000) T2(w,c0,v0,a0) Fc1(v0,c0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_v0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == 0 .and. & 
IEOR(s_w,s_c0) == IEOR(s_v0,s_a0) .and. &
IEOR(s_v0,s_c0) == 0) then

if(psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- T2(w,c0,v0,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_c0, i_v0) =  &
  T2_(s_v0, s_c0, s_w)%array(i_v0, i_c0, i_w)
end do
end do
end do
! Z2 <-- Fc1(v0,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_v0) =  &
  Fc1_(s_c0, s_v0)%array(i_c0, i_v0)
end do
end do

! Z3 <-- W0(w,a0) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w),&
                     1,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w))

! W0(w,a0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
W0_(s_w)%array(i_w) = &
    W0_(s_w)%array(i_w) &
  + Z3_(i_w)
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w) * &
                1 * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no0_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no0_x1_type0_noeri &
  (sa0, ia0, sj, ij, W0, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0, sj, ij
real(kind=8), intent(inout) :: W0(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa0

call set_symblock_Xc(sleft, W0, nir, nsym, psym) ! -> Xc (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no0_x1_type0_noeri &
  (sa0, ia0, sj, ij, Xc, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no0_x1_type0_noeri



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
subroutine g_sigma_cooo_ccov_no0_x1_type0_noeri &
  (s_a0, i_a0, s_j, i_j, W0_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W0_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (    1.00000000) D2(k,j,a0,i) W0(w,a0) 
do s_k = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == IEOR(s_a0,s_i) .and. &
IEOR(s_w,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0) then

! Z1 <-- D2(k,j,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i) =  &
  D2_(s_i, s_a0, s_j, s_k)%array(i_i, i_a0, i_j, i_k)
end do
end do
! Z2 <-- W0(w,a0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z2_(i_w) =  &
  W0_(s_w)%array(i_w)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     1,&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no0_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no1_x0_type0_noeri &
  (sv0, iv0, T2, W1, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), W1(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W1, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no1_x0_type0_noeri &
  (sv0, iv0, av2_i, Xca, fc1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no1_x0_type0_noeri



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
subroutine g_sigma_cooo_ccov_no1_x0_type0_noeri &
  (s_v0, i_v0, T2_, W1_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W1_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_a0, i_a0
! W1(w,a0) += (    1.00000000) T2(w,c0,a0,v0) Fc1(v0,c0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == 0 .and. & 
IEOR(s_w,s_c0) == IEOR(s_a0,s_v0) .and. &
IEOR(s_v0,s_c0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(w,c0,a0,v0) 
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
! Z2 <-- Fc1(v0,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0) =  &
  Fc1_(s_c0, s_v0)%array(i_c0, i_v0)
end do

! Z3 <-- W1(w,a0) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0))

! W1(w,a0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W1_(s_a0, s_w)%array(i_a0, i_w) = &
    W1_(s_a0, s_w)%array(i_a0, i_w) &
  + Z3_(i_w, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no1_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no1_x1_type0_noeri &
  (sj, ij, W1, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W1(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W1, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no1_x1_type0_noeri &
  (sj, ij, Xca, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no1_x1_type0_noeri



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
subroutine g_sigma_cooo_ccov_no1_x1_type0_noeri &
  (s_j, i_j, W1_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W1_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a0, i_a0, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (   -2.00000000) D2(k,j,a0,i) W1(w,a0) 
do s_k = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == IEOR(s_a0,s_i) .and. &
IEOR(s_w,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(k,j,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a0) =  &
  D2_(s_i, s_a0, s_j, s_k)%array(i_i, i_a0, i_j, i_k)
end do
end do
end do
! Z2 <-- W1(w,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  W1_(s_a0, s_w)%array(i_a0, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no1_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no2_x0_type0_noeri &
  (sv0, iv0, T2, W2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), W2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W2, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no2_x0_type0_noeri &
  (sv0, iv0, av2_i, Xca, fc1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no2_x0_type0_noeri



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
subroutine g_sigma_cooo_ccov_no2_x0_type0_noeri &
  (s_v0, i_v0, T2_, W2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W2_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_i, i_i
! W2(w,i) += (    1.00000000) T2(w,c0,i,v0) Fc1(v0,c0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == 0 .and. & 
IEOR(s_w,s_c0) == IEOR(s_i,s_v0) .and. &
IEOR(s_v0,s_c0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(w,c0,i,v0) 
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
! Z2 <-- Fc1(v0,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0) =  &
  Fc1_(s_c0, s_v0)%array(i_c0, i_v0)
end do

! Z3 <-- W2(w,i) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! W2(w,i)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W2_(s_i, s_w)%array(i_i, i_w) = &
    W2_(s_i, s_w)%array(i_i, i_w) &
  + Z3_(i_w, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                1 * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no2_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no2_x1_type0_noeri &
  (sj, ij, W2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W2, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no2_x1_type0_noeri &
  (sj, ij, Xca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no2_x1_type0_noeri



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
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_cooo_ccov_no2_x1_type0_noeri &
  (s_j, i_j, W2_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W2_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_w, i_w, s_i, i_i
! S2(w,k,i,j) += (    4.00000000) D1(k,j) W2(w,i) 
do s_k = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == 0 .and. &
IEOR(s_w,s_i) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0) then

! Z1 <-- W2(w,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i) =  &
  W2_(s_i, s_w)%array(i_i, i_w)
end do
end do
! Z2 <-- D1(k,j) 
allocate(Z2_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_k) =  &
  D1_(s_j, s_k)%array(i_j, i_k)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     1,&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no2_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no3_x0_type0_noeri &
  (si, ii, T2, W3, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: si, ii
real(kind=8), intent(inout) :: T2(*), W3(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(si, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = si

call set_symblock_Xc(sleft, W3, nir, nsym, psym) ! -> Xc (allocate) 
call g_sigma_cooo_ccov_no3_x0_type0_noeri &
  (si, ii, av2_i, Xc, fc1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no3_x0_type0_noeri



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
subroutine g_sigma_cooo_ccov_no3_x0_type0_noeri &
  (s_i, i_i, T2_, W3_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_i, s_i
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W3_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_v0, i_v0
! W3(w,i) += (    1.00000000) T2(w,c0,v0,i) Fc1(v0,c0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_v0 = 0, nir-1
if( &
IEOR(s_w,s_i) == 0 .and. & 
IEOR(s_w,s_c0) == IEOR(s_v0,s_i) .and. &
IEOR(s_v0,s_c0) == 0) then

if(psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- T2(w,c0,v0,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_c0, i_v0) =  &
  T2_(s_v0, s_c0, s_w)%array(i_v0, i_c0, i_w)
end do
end do
end do
! Z2 <-- Fc1(v0,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_v0) =  &
  Fc1_(s_c0, s_v0)%array(i_c0, i_v0)
end do
end do

! Z3 <-- W3(w,i) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w),&
                     1,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w))

! W3(w,i)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
W3_(s_w)%array(i_w) = &
    W3_(s_w)%array(i_w) &
  + Z3_(i_w)
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w) * &
                1 * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no3_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no3_x1_type0_noeri &
  (si, ii, sj, ij, W3, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: si, ii, sj, ij
real(kind=8), intent(inout) :: W3(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = si

call set_symblock_Xc(sleft, W3, nir, nsym, psym) ! -> Xc (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no3_x1_type0_noeri &
  (si, ii, sj, ij, Xc, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no3_x1_type0_noeri



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
subroutine g_sigma_cooo_ccov_no3_x1_type0_noeri &
  (s_i, i_i, s_j, i_j, W3_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
integer, intent(in) :: i_i, s_i
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W3_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_w, i_w
! S2(w,k,i,j) += (   -2.00000000) D1(k,j) W3(w,i) 
do s_k = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == 0 .and. &
IEOR(s_w,s_i) == 0) then

if(psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0) then

! Z1 <-- D1(k,j) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k) =  &
  D1_(s_j, s_k)%array(i_j, i_k)
end do
! Z2 <-- W3(w,i) 
allocate(Z2_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z2_(i_w) =  &
  W3_(s_w)%array(i_w)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_w),&
                     1,&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_w) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no3_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no4_x0_type0_noeri &
  (sj, ij, T2, W4, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: T2(*), W4(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sj, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sj

call set_symblock_Xc(sleft, W4, nir, nsym, psym) ! -> Xc (allocate) 
call g_sigma_cooo_ccov_no4_x0_type0_noeri &
  (sj, ij, av2_i, Xc, fc1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no4_x0_type0_noeri



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
subroutine g_sigma_cooo_ccov_no4_x0_type0_noeri &
  (s_j, i_j, T2_, W4_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W4_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_w, i_w, s_v0, i_v0
! W4(w,j) += (    1.00000000) T2(c0,w,v0,j) Fc1(v0,c0) 
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_v0 = 0, nir-1
if( &
IEOR(s_w,s_j) == 0 .and. & 
IEOR(s_c0,s_w) == IEOR(s_v0,s_j) .and. &
IEOR(s_v0,s_c0) == 0) then

if(psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- T2(c0,w,v0,j) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_c0, i_v0) =  &
  T2_(s_v0, s_w, s_c0)%array(i_v0, i_w, i_c0)
end do
end do
end do
! Z2 <-- Fc1(v0,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_v0) =  &
  Fc1_(s_c0, s_v0)%array(i_c0, i_v0)
end do
end do

! Z3 <-- W4(w,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w),&
                     1,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w))

! W4(w,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
W4_(s_w)%array(i_w) = &
    W4_(s_w)%array(i_w) &
  + Z3_(i_w)
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w) * &
                1 * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no4_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no4_x1_type0_noeri &
  (sj, ij, W4, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W4(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sj

call set_symblock_Xc(sleft, W4, nir, nsym, psym) ! -> Xc (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no4_x1_type0_noeri &
  (sj, ij, Xc, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no4_x1_type0_noeri



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
subroutine g_sigma_cooo_ccov_no4_x1_type0_noeri &
  (s_j, i_j, W4_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W4_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (   -2.00000000) D1(k,i) W4(w,j) 
do s_k = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_i) == 0 .and. &
IEOR(s_w,s_j) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0) then

! Z1 <-- D1(k,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i) =  &
  D1_(s_i, s_k)%array(i_i, i_k)
end do
end do
! Z2 <-- W4(w,j) 
allocate(Z2_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z2_(i_w) =  &
  W4_(s_w)%array(i_w)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     1,&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no4_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no5_x0_type0_noeri &
  (sj, ij, T2, W5, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: T2(*), W5(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sj, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sj

call set_symblock_Xc(sleft, W5, nir, nsym, psym) ! -> Xc (allocate) 
call g_sigma_cooo_ccov_no5_x0_type0_noeri &
  (sj, ij, av2_i, Xc, fc1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no5_x0_type0_noeri



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
subroutine g_sigma_cooo_ccov_no5_x0_type0_noeri &
  (s_j, i_j, T2_, W5_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W5_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_v0, i_v0
! W5(w,j) += (    1.00000000) T2(w,c0,v0,j) Fc1(v0,c0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_v0 = 0, nir-1
if( &
IEOR(s_w,s_j) == 0 .and. & 
IEOR(s_w,s_c0) == IEOR(s_v0,s_j) .and. &
IEOR(s_v0,s_c0) == 0) then

if(psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- T2(w,c0,v0,j) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_c0, i_v0) =  &
  T2_(s_v0, s_c0, s_w)%array(i_v0, i_c0, i_w)
end do
end do
end do
! Z2 <-- Fc1(v0,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_v0) =  &
  Fc1_(s_c0, s_v0)%array(i_c0, i_v0)
end do
end do

! Z3 <-- W5(w,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w),&
                     1,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w))

! W5(w,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
W5_(s_w)%array(i_w) = &
    W5_(s_w)%array(i_w) &
  + Z3_(i_w)
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w) * &
                1 * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no5_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no5_x1_type0_noeri &
  (sj, ij, W5, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W5(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sj

call set_symblock_Xc(sleft, W5, nir, nsym, psym) ! -> Xc (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no5_x1_type0_noeri &
  (sj, ij, Xc, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no5_x1_type0_noeri



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
subroutine g_sigma_cooo_ccov_no5_x1_type0_noeri &
  (s_j, i_j, W5_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W5_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (    1.00000000) D1(k,i) W5(w,j) 
do s_k = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_i) == 0 .and. &
IEOR(s_w,s_j) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0) then

! Z1 <-- D1(k,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i) =  &
  D1_(s_i, s_k)%array(i_i, i_k)
end do
end do
! Z2 <-- W5(w,j) 
allocate(Z2_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z2_(i_w) =  &
  W5_(s_w)%array(i_w)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     1,&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no5_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no0_x0_type1_eri_c &
  (sj, ij, sw, iw, T2, V2, W24, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sw, iw
real(kind=8), intent(inout) :: T2(*), V2(*), W24
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sj, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sw,sj)

call g_sigma_cooo_ccov_no0_x0_type1_eri_c &
  (sj, ij, sw, iw, av2_i, h2_i, W24, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_cooo_ccov_no0_x0_type1_eri_c



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
subroutine g_sigma_cooo_ccov_no0_x0_type1_eri_c &
  (s_j, i_j, s_w, i_w, T2_, V2_, W24_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W24_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_v0, i_v0, s_c0, i_c0
! W24(w,j) += (    1.00000000) V2(w,c1,v0,c0) T2(c0,c1,v0,j) 
do s_c1 = 0, nir-1
do s_v0 = 0, nir-1
do s_c0 = 0, nir-1
if( &
IEOR(s_w,s_j) == 0 .and. & 
IEOR(s_w,s_c1) == IEOR(s_v0,s_c0) .and. &
IEOR(s_c0,s_c1) == IEOR(s_v0,s_j)) then

if(psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(w,c1,v0,c0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z1_(i_c1, i_v0, i_c0) =  &
  V2_(s_c0, s_v0, s_c1)%array(i_c0, i_v0, i_c1)
end do
end do
end do
! Z2 <-- T2(c0,c1,v0,j) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z2_(i_c1, i_v0, i_c0) =  &
  T2_(s_v0, s_c1, s_c0)%array(i_v0, i_c1, i_c0)
end do
end do
end do

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', 1,&
                     1,&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     1,&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     1)

! W24(w,j)  <-- Z3
W24_ = &
    W24_ &
  + Z3_

! Flop count
flops = flops + 1 * &
                1 * &
                psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no0_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no0_x1_type1_eri_c &
  (sj, ij, sw, iw, W24, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sw, iw
real(kind=8), intent(inout) :: W24, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sw,sj)

call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no0_x1_type1_eri_c &
  (sj, ij, sw, iw, W24, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)

end subroutine g_if_sigma_cooo_ccov_no0_x1_type1_eri_c



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
subroutine g_sigma_cooo_ccov_no0_x1_type1_eri_c &
  (s_j, i_j, s_w, i_w, W24_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W24_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_i, i_i
! S2(w,k,i,j) += (    2.00000000) D1(k,i) W24(w,j) 
do s_k = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_i) == 0 .and. &
IEOR(s_w,s_j) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0) then

! Z1 <-- D1(k,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i) =  &
  D1_(s_i, s_k)%array(i_i, i_k)
end do
end do
! Z2 <-- W24(w,j) 
Z2_ =  &
  W24_

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     1,&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                1 * &
                1 * 2.0d+00

deallocate(Z1_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no0_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no1_x0_type1_eri_c &
  (sj, ij, sw, iw, T2, V2, W25, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sw, iw
real(kind=8), intent(inout) :: T2(*), V2(*), W25
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sj, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sw,sj)

call g_sigma_cooo_ccov_no1_x0_type1_eri_c &
  (sj, ij, sw, iw, av2_i, h2_i, W25, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_cooo_ccov_no1_x0_type1_eri_c



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
subroutine g_sigma_cooo_ccov_no1_x0_type1_eri_c &
  (s_j, i_j, s_w, i_w, T2_, V2_, W25_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W25_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_v0, i_v0, s_c0, i_c0
! W25(w,j) += (    1.00000000) V2(w,c1,v0,c0) T2(c1,c0,v0,j) 
do s_c1 = 0, nir-1
do s_v0 = 0, nir-1
do s_c0 = 0, nir-1
if( &
IEOR(s_w,s_j) == 0 .and. & 
IEOR(s_w,s_c1) == IEOR(s_v0,s_c0) .and. &
IEOR(s_c1,s_c0) == IEOR(s_v0,s_j)) then

if(psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(w,c1,v0,c0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z1_(i_c1, i_v0, i_c0) =  &
  V2_(s_c0, s_v0, s_c1)%array(i_c0, i_v0, i_c1)
end do
end do
end do
! Z2 <-- T2(c1,c0,v0,j) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z2_(i_c1, i_v0, i_c0) =  &
  T2_(s_v0, s_c0, s_c1)%array(i_v0, i_c0, i_c1)
end do
end do
end do

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', 1,&
                     1,&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     1,&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     1)

! W25(w,j)  <-- Z3
W25_ = &
    W25_ &
  + Z3_

! Flop count
flops = flops + 1 * &
                1 * &
                psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no1_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no1_x1_type1_eri_c &
  (sj, ij, sw, iw, W25, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sw, iw
real(kind=8), intent(inout) :: W25, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sw,sj)

call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no1_x1_type1_eri_c &
  (sj, ij, sw, iw, W25, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)

end subroutine g_if_sigma_cooo_ccov_no1_x1_type1_eri_c



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
subroutine g_sigma_cooo_ccov_no1_x1_type1_eri_c &
  (s_j, i_j, s_w, i_w, W25_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W25_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_i, i_i
! S2(w,k,i,j) += (   -1.00000000) D1(k,i) W25(w,j) 
do s_k = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_i) == 0 .and. &
IEOR(s_w,s_j) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0) then

! Z1 <-- D1(k,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i) =  &
  D1_(s_i, s_k)%array(i_i, i_k)
end do
end do
! Z2 <-- W25(w,j) 
Z2_ =  &
  W25_

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     1,&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                1 * &
                1 * 2.0d+00

deallocate(Z1_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no1_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no0_x0_type1_eri_o &
  (sa0, ia0, sj, ij, T2, V2, W22, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0, sj, ij
real(kind=8), intent(inout) :: T2(*), V2(*), W22(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sj, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sj)

call set_symblock_Xca(sleft, W22, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no0_x0_type1_eri_o &
  (sa0, ia0, sj, ij, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no0_x0_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no0_x0_type1_eri_o &
  (s_a0, i_a0, s_j, i_j, T2_, V2_, W22_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
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

integer :: s_a1, i_a1, s_v0, i_v0, s_c0, i_c0, s_w, i_w
! W22(w,a0,j,a1) += (    1.00000000) V2(j,a1,v0,c0) T2(w,c0,v0,a0) 
do s_a1 = 0, nir-1
do s_v0 = 0, nir-1
do s_c0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_j,s_a1) .and. & 
IEOR(s_j,s_a1) == IEOR(s_v0,s_c0) .and. &
IEOR(s_w,s_c0) == IEOR(s_v0,s_a0)) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(j,a1,v0,c0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_v0, i_c0) =  &
  V2_(s_c0, s_v0, s_a1)%array(i_c0, i_v0, i_a1)
end do
end do
end do
! Z2 <-- T2(w,c0,v0,a0) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_c0, i_w) =  &
  T2_(s_v0, s_c0, s_w)%array(i_v0, i_c0, i_w)
end do
end do
end do

! Z3 <-- W22(w,a0,j,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W22(w,a0,j,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W22_(s_a1, s_w)%array(i_a1, i_w) = &
    W22_(s_a1, s_w)%array(i_a1, i_w) &
  + Z3_(i_a1, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no0_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no0_x1_type1_eri_o &
  (sa0, ia0, sj, ij, W22, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0, sj, ij
real(kind=8), intent(inout) :: W22(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa0,sj)

call set_symblock_Xca(sleft, W22, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no0_x1_type1_eri_o &
  (sa0, ia0, sj, ij, Xca, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no0_x1_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no0_x1_type1_eri_o &
  (s_a0, i_a0, s_j, i_j, W22_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W22_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a1, i_a1, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (    1.00000000) D2(k,a1,a0,i) W22(w,a0,j,a1) 
do s_k = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_a1) == IEOR(s_a0,s_i) .and. &
IEOR(s_w,s_a0) == IEOR(s_j,s_a1)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(k,a1,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a1) =  &
  D2_(s_i, s_a0, s_a1, s_k)%array(i_i, i_a0, i_a1, i_k)
end do
end do
end do
! Z2 <-- W22(w,a0,j,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w) =  &
  W22_(s_a1, s_w)%array(i_a1, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no0_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no1_x0_type1_eri_o &
  (sj, ij, sv0, iv0, T2, V2, W23, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W23(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sj, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sj

call set_symblock_Xcaa(sleft, W23, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_cooo_ccov_no1_x0_type1_eri_o &
  (sj, ij, sv0, iv0, av2_i, h2_i, Xcaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no1_x0_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no1_x0_type1_eri_o &
  (s_j, i_j, s_v0, i_v0, T2_, V2_, W23_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W23_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_c0, i_c0, s_w, i_w, s_a0, i_a0
! W23(w,a0,j,a1) += (    1.00000000) V2(j,a1,v0,c0) T2(w,c0,a0,v0) 
do s_a1 = 0, nir-1
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_j,s_a1) .and. & 
IEOR(s_j,s_a1) == IEOR(s_v0,s_c0) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(j,a1,v0,c0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_c0) =  &
  V2_(s_c0, s_v0, s_a1)%array(i_c0, i_v0, i_a1)
end do
end do
! Z2 <-- T2(w,c0,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w, i_a0) =  &
  T2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- W23(w,a0,j,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W23(w,a0,j,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W23_(s_a1, s_a0, s_w)%array(i_a1, i_a0, i_w) = &
    W23_(s_a1, s_a0, s_w)%array(i_a1, i_a0, i_w) &
  + Z3_(i_a1, i_w, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no1_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no1_x1_type1_eri_o &
  (sj, ij, W23, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W23(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sj

call set_symblock_Xcaa(sleft, W23, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no1_x1_type1_eri_o &
  (sj, ij, Xcaa, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no1_x1_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no1_x1_type1_eri_o &
  (s_j, i_j, W23_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W23_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a1, i_a1, s_a0, i_a0, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (   -2.00000000) D2(k,a1,a0,i) W23(w,a0,j,a1) 
do s_k = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_a1) == IEOR(s_a0,s_i) .and. &
IEOR(s_w,s_a0) == IEOR(s_j,s_a1)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(k,a1,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a1, i_a0) =  &
  D2_(s_i, s_a0, s_a1, s_k)%array(i_i, i_a0, i_a1, i_k)
end do
end do
end do
end do
! Z2 <-- W23(w,a0,j,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_w) =  &
  W23_(s_a1, s_a0, s_w)%array(i_a1, i_a0, i_w)
end do
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no1_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no2_x0_type1_eri_o &
  (sa0, ia0, sj, ij, T2, V2, W30, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0, sj, ij
real(kind=8), intent(inout) :: T2(*), V2(*), W30(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sj, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sj)

call set_symblock_Xca(sleft, W30, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no2_x0_type1_eri_o &
  (sa0, ia0, sj, ij, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no2_x0_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no2_x0_type1_eri_o &
  (s_a0, i_a0, s_j, i_j, T2_, V2_, W30_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W30_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_c0, i_c0, s_a1, i_a1, s_w, i_w
! W30(w,a0,j,a1) += (    1.00000000) V2(j,v0,c0,a1) T2(w,c0,v0,a0) 
do s_v0 = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_j,s_a1) .and. & 
IEOR(s_j,s_v0) == IEOR(s_c0,s_a1) .and. &
IEOR(s_w,s_c0) == IEOR(s_v0,s_a0)) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(j,v0,c0,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_v0, i_c0) =  &
  V2_(s_a1, s_c0, s_v0)%array(i_a1, i_c0, i_v0)
end do
end do
end do
! Z2 <-- T2(w,c0,v0,a0) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_c0, i_w) =  &
  T2_(s_v0, s_c0, s_w)%array(i_v0, i_c0, i_w)
end do
end do
end do

! Z3 <-- W30(w,a0,j,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W30(w,a0,j,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W30_(s_a1, s_w)%array(i_a1, i_w) = &
    W30_(s_a1, s_w)%array(i_a1, i_w) &
  + Z3_(i_a1, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no2_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no2_x1_type1_eri_o &
  (sa0, ia0, sj, ij, W30, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0, sj, ij
real(kind=8), intent(inout) :: W30(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa0,sj)

call set_symblock_Xca(sleft, W30, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no2_x1_type1_eri_o &
  (sa0, ia0, sj, ij, Xca, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no2_x1_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no2_x1_type1_eri_o &
  (s_a0, i_a0, s_j, i_j, W30_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a0, s_a0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W30_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_i, i_i, s_a1, i_a1, s_w, i_w
! S2(w,k,i,j) += (    1.00000000) D2(k,i,a0,a1) W30(w,a0,j,a1) 
do s_k = 0, nir-1
do s_i = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_i) == IEOR(s_a0,s_a1) .and. &
IEOR(s_w,s_a0) == IEOR(s_j,s_a1)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(k,i,a0,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a1) =  &
  D2_(s_a1, s_a0, s_i, s_k)%array(i_a1, i_a0, i_i, i_k)
end do
end do
end do
! Z2 <-- W30(w,a0,j,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w) =  &
  W30_(s_a1, s_w)%array(i_a1, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no2_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no3_x0_type1_eri_o &
  (sj, ij, sv0, iv0, T2, V2, W31, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W31(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sj, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sj

call set_symblock_Xcaa(sleft, W31, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_cooo_ccov_no3_x0_type1_eri_o &
  (sj, ij, sv0, iv0, av2_i, h2_i, Xcaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no3_x0_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no3_x0_type1_eri_o &
  (s_j, i_j, s_v0, i_v0, T2_, V2_, W31_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
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

integer :: s_c0, i_c0, s_a1, i_a1, s_w, i_w, s_a0, i_a0
! W31(w,a0,j,a1) += (    1.00000000) V2(j,v0,c0,a1) T2(w,c0,a0,v0) 
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_j,s_a1) .and. & 
IEOR(s_j,s_v0) == IEOR(s_c0,s_a1) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(j,v0,c0,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_c0) =  &
  V2_(s_a1, s_c0, s_v0)%array(i_a1, i_c0, i_v0)
end do
end do
! Z2 <-- T2(w,c0,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w, i_a0) =  &
  T2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- W31(w,a0,j,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W31(w,a0,j,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W31_(s_a1, s_a0, s_w)%array(i_a1, i_a0, i_w) = &
    W31_(s_a1, s_a0, s_w)%array(i_a1, i_a0, i_w) &
  + Z3_(i_a1, i_w, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no3_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no3_x1_type1_eri_o &
  (sj, ij, W31, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W31(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sj

call set_symblock_Xcaa(sleft, W31, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no3_x1_type1_eri_o &
  (sj, ij, Xcaa, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no3_x1_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no3_x1_type1_eri_o &
  (s_j, i_j, W31_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W31_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a1, i_a1, s_a0, i_a0, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (    1.00000000) D2(k,a1,a0,i) W31(w,a0,j,a1) 
do s_k = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_a1) == IEOR(s_a0,s_i) .and. &
IEOR(s_w,s_a0) == IEOR(s_j,s_a1)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(k,a1,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a1, i_a0) =  &
  D2_(s_i, s_a0, s_a1, s_k)%array(i_i, i_a0, i_a1, i_k)
end do
end do
end do
end do
! Z2 <-- W31(w,a0,j,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_w) =  &
  W31_(s_a1, s_a0, s_w)%array(i_a1, i_a0, i_w)
end do
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no3_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no4_x0_type1_eri_o &
  (sj, ij, sv0, iv0, T2, V2, W48, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W48(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sj, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sj

call set_symblock_Xc(sleft, W48, nir, nsym, psym) ! -> Xc (allocate) 
call g_sigma_cooo_ccov_no4_x0_type1_eri_o &
  (sj, ij, sv0, iv0, av2_i, h2_i, Xc, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no4_x0_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no4_x0_type1_eri_o &
  (s_j, i_j, s_v0, i_v0, T2_, V2_, W48_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W48_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_c0, i_c0, s_w, i_w
! W48(w,j) += (    1.00000000) V2(j,a0,v0,c0) T2(w,c0,a0,v0) 
do s_a0 = 0, nir-1
do s_c0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_j) == 0 .and. & 
IEOR(s_j,s_a0) == IEOR(s_v0,s_c0) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(w,c0,a0,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_c0, i_a0) =  &
  T2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do
! Z2 <-- V2(j,a0,v0,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0) =  &
  V2_(s_c0, s_v0, s_a0)%array(i_c0, i_v0, i_a0)
end do
end do

! Z3 <-- W48(w,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w),&
                     1,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w))

! W48(w,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
W48_(s_w)%array(i_w) = &
    W48_(s_w)%array(i_w) &
  + Z3_(i_w)
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w) * &
                1 * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no4_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no4_x1_type1_eri_o &
  (sj, ij, W48, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W48(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sj

call set_symblock_Xc(sleft, W48, nir, nsym, psym) ! -> Xc (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no4_x1_type1_eri_o &
  (sj, ij, Xc, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no4_x1_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no4_x1_type1_eri_o &
  (s_j, i_j, W48_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W48_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (   -2.00000000) D1(k,i) W48(w,j) 
do s_k = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_i) == 0 .and. &
IEOR(s_w,s_j) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0) then

! Z1 <-- D1(k,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i) =  &
  D1_(s_i, s_k)%array(i_i, i_k)
end do
end do
! Z2 <-- W48(w,j) 
allocate(Z2_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z2_(i_w) =  &
  W48_(s_w)%array(i_w)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     1,&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no4_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no5_x0_type1_eri_o &
  (sa0, ia0, sj, ij, T2, V2, W49, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0, sj, ij
real(kind=8), intent(inout) :: T2(*), V2(*), W49(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sj, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sj

call set_symblock_Xc(sleft, W49, nir, nsym, psym) ! -> Xc (allocate) 
call g_sigma_cooo_ccov_no5_x0_type1_eri_o &
  (sa0, ia0, sj, ij, av2_i, h2_i, Xc, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no5_x0_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no5_x0_type1_eri_o &
  (s_a0, i_a0, s_j, i_j, T2_, V2_, W49_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W49_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_c0, i_c0, s_w, i_w
! W49(w,j) += (    1.00000000) V2(j,a0,v0,c0) T2(w,c0,v0,a0) 
do s_v0 = 0, nir-1
do s_c0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_j) == 0 .and. & 
IEOR(s_j,s_a0) == IEOR(s_v0,s_c0) .and. &
IEOR(s_w,s_c0) == IEOR(s_v0,s_a0)) then

if(psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- T2(w,c0,v0,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_c0, i_v0) =  &
  T2_(s_v0, s_c0, s_w)%array(i_v0, i_c0, i_w)
end do
end do
end do
! Z2 <-- V2(j,a0,v0,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_v0) =  &
  V2_(s_c0, s_v0, s_a0)%array(i_c0, i_v0, i_a0)
end do
end do

! Z3 <-- W49(w,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w),&
                     1,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w))

! W49(w,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
W49_(s_w)%array(i_w) = &
    W49_(s_w)%array(i_w) &
  + Z3_(i_w)
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w) * &
                1 * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no5_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no5_x1_type1_eri_o &
  (sj, ij, W49, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W49(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sj

call set_symblock_Xc(sleft, W49, nir, nsym, psym) ! -> Xc (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no5_x1_type1_eri_o &
  (sj, ij, Xc, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no5_x1_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no5_x1_type1_eri_o &
  (s_j, i_j, W49_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W49_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (    1.00000000) D1(k,i) W49(w,j) 
do s_k = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_i) == 0 .and. &
IEOR(s_w,s_j) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0) then

! Z1 <-- D1(k,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i) =  &
  D1_(s_i, s_k)%array(i_i, i_k)
end do
end do
! Z2 <-- W49(w,j) 
allocate(Z2_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z2_(i_w) =  &
  W49_(s_w)%array(i_w)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     1,&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no5_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no6_x0_type1_eri_o &
  (sa0, ia0, sj, ij, T2, V2, W50, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0, sj, ij
real(kind=8), intent(inout) :: T2(*), V2(*), W50(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sj, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sj

call set_symblock_Xc(sleft, W50, nir, nsym, psym) ! -> Xc (allocate) 
call g_sigma_cooo_ccov_no6_x0_type1_eri_o &
  (sa0, ia0, sj, ij, av2_i, h2_i, Xc, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no6_x0_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no6_x0_type1_eri_o &
  (s_a0, i_a0, s_j, i_j, T2_, V2_, W50_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W50_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_c0, i_c0, s_w, i_w
! W50(w,j) += (    1.00000000) V2(j,v0,c0,a0) T2(w,c0,v0,a0) 
do s_v0 = 0, nir-1
do s_c0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_j) == 0 .and. & 
IEOR(s_j,s_v0) == IEOR(s_c0,s_a0) .and. &
IEOR(s_w,s_c0) == IEOR(s_v0,s_a0)) then

if(psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- T2(w,c0,v0,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_c0, i_v0) =  &
  T2_(s_v0, s_c0, s_w)%array(i_v0, i_c0, i_w)
end do
end do
end do
! Z2 <-- V2(j,v0,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_v0) =  &
  V2_(s_a0, s_c0, s_v0)%array(i_a0, i_c0, i_v0)
end do
end do

! Z3 <-- W50(w,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w),&
                     1,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w))

! W50(w,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
W50_(s_w)%array(i_w) = &
    W50_(s_w)%array(i_w) &
  + Z3_(i_w)
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w) * &
                1 * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no6_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no6_x1_type1_eri_o &
  (sj, ij, W50, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W50(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sj

call set_symblock_Xc(sleft, W50, nir, nsym, psym) ! -> Xc (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no6_x1_type1_eri_o &
  (sj, ij, Xc, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no6_x1_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no6_x1_type1_eri_o &
  (s_j, i_j, W50_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W50_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (   -2.00000000) D1(k,i) W50(w,j) 
do s_k = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_i) == 0 .and. &
IEOR(s_w,s_j) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0) then

! Z1 <-- D1(k,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i) =  &
  D1_(s_i, s_k)%array(i_i, i_k)
end do
end do
! Z2 <-- W50(w,j) 
allocate(Z2_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z2_(i_w) =  &
  W50_(s_w)%array(i_w)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     1,&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no6_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no7_x0_type1_eri_o &
  (sj, ij, sv0, iv0, T2, V2, W51, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W51(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sj, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sj

call set_symblock_Xc(sleft, W51, nir, nsym, psym) ! -> Xc (allocate) 
call g_sigma_cooo_ccov_no7_x0_type1_eri_o &
  (sj, ij, sv0, iv0, av2_i, h2_i, Xc, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no7_x0_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no7_x0_type1_eri_o &
  (s_j, i_j, s_v0, i_v0, T2_, V2_, W51_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W51_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a0, i_a0, s_w, i_w
! W51(w,j) += (    1.00000000) V2(j,v0,c0,a0) T2(w,c0,a0,v0) 
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_j) == 0 .and. & 
IEOR(s_j,s_v0) == IEOR(s_c0,s_a0) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(w,c0,a0,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_c0, i_a0) =  &
  T2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do
! Z2 <-- V2(j,v0,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0) =  &
  V2_(s_a0, s_c0, s_v0)%array(i_a0, i_c0, i_v0)
end do
end do

! Z3 <-- W51(w,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w),&
                     1,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w))

! W51(w,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
W51_(s_w)%array(i_w) = &
    W51_(s_w)%array(i_w) &
  + Z3_(i_w)
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w) * &
                1 * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no7_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no7_x1_type1_eri_o &
  (sj, ij, W51, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W51(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sj

call set_symblock_Xc(sleft, W51, nir, nsym, psym) ! -> Xc (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no7_x1_type1_eri_o &
  (sj, ij, Xc, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xc)

end subroutine g_if_sigma_cooo_ccov_no7_x1_type1_eri_o



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
subroutine g_sigma_cooo_ccov_no7_x1_type1_eri_o &
  (s_j, i_j, W51_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock1), intent(inout) :: W51_(0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (    1.00000000) D1(k,i) W51(w,j) 
do s_k = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_i) == 0 .and. &
IEOR(s_w,s_j) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0) then

! Z1 <-- D1(k,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i) =  &
  D1_(s_i, s_k)%array(i_i, i_k)
end do
end do
! Z2 <-- W51(w,j) 
allocate(Z2_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z2_(i_w) =  &
  W51_(s_w)%array(i_w)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     1,&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no7_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no0_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W6, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W6(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W6, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no0_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no0_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no0_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W6_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W6_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_w, i_w, s_c1, i_c1, s_a0, i_a0
! W6(w,a0) += (    1.00000000) V2(v0,c0,w,c1) T2(c0,c1,a0,v0) 
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_c1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == 0 .and. & 
IEOR(s_v0,s_c0) == IEOR(s_w,s_c1) .and. &
IEOR(s_c0,s_c1) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1) > 0) then

! Z1 <-- T2(c0,c1,a0,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1)))
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_c0, i_c1) =  &
  T2_(s_a0, s_c1, s_c0)%array(i_a0, i_c1, i_c0)
end do
end do
end do
! Z2 <-- V2(v0,c0,w,c1) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_c1, i_w) =  &
  V2_(s_c1, s_w, s_c0)%array(i_c1, i_w, i_c0)
end do
end do
end do

! Z3 <-- W6(w,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W6(w,a0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W6_(s_a0, s_w)%array(i_a0, i_w) = &
    W6_(s_a0, s_w)%array(i_a0, i_w) &
  + Z3_(i_a0, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no0_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no1_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W7, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W7(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W7, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no1_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no1_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no1_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W7_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W7_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_w, i_w, s_c1, i_c1, s_a0, i_a0
! W7(w,a0) += (    1.00000000) V2(v0,c0,w,c1) T2(c1,c0,a0,v0) 
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_c1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == 0 .and. & 
IEOR(s_v0,s_c0) == IEOR(s_w,s_c1) .and. &
IEOR(s_c1,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c1,c0,a0,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z1_(i_a0, i_c1, i_c0) =  &
  T2_(s_a0, s_c0, s_c1)%array(i_a0, i_c0, i_c1)
end do
end do
end do
! Z2 <-- V2(v0,c0,w,c1) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z2_(i_c1, i_c0, i_w) =  &
  V2_(s_c1, s_w, s_c0)%array(i_c1, i_w, i_c0)
end do
end do
end do

! Z3 <-- W7(w,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W7(w,a0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W7_(s_a0, s_w)%array(i_a0, i_w) = &
    W7_(s_a0, s_w)%array(i_a0, i_w) &
  + Z3_(i_a0, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no1_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no2_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W8, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W8(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xcaaa(sleft, W8, nir, nsym, psym) ! -> Xcaaa (allocate) 
call g_sigma_cooo_ccov_no2_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xcaaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no2_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no2_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W8_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W8_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_c0, i_c0, s_a1, i_a1, s_w, i_w, s_a0, i_a0
! W8(w,a0,a2,a1) += (    1.00000000) V2(v0,a2,c0,a1) T2(c0,w,a0,v0) 
do s_a2 = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_a2,s_a1) .and. & 
IEOR(s_v0,s_a2) == IEOR(s_c0,s_a1) .and. &
IEOR(s_c0,s_w) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(v0,a2,c0,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_a1, i_c0) =  &
  V2_(s_a1, s_c0, s_a2)%array(i_a1, i_c0, i_a2)
end do
end do
end do
! Z2 <-- T2(c0,w,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w, i_a0) =  &
  T2_(s_a0, s_w, s_c0)%array(i_a0, i_w, i_c0)
end do
end do
end do

! Z3 <-- W8(w,a0,a2,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1))

! W8(w,a0,a2,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W8_(s_a1, s_a2, s_a0, s_w)%array(i_a1, i_a2, i_a0, i_w) = &
    W8_(s_a1, s_a2, s_a0, s_w)%array(i_a1, i_a2, i_a0, i_w) &
  + Z3_(i_a2, i_a1, i_w, i_a0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no2_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no3_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W9, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W9(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xcaaa(sleft, W9, nir, nsym, psym) ! -> Xcaaa (allocate) 
call g_sigma_cooo_ccov_no3_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xcaaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no3_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no3_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W9_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W9_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_c0, i_c0, s_a1, i_a1, s_w, i_w, s_a0, i_a0
! W9(w,a0,a2,a1) += (    1.00000000) V2(v0,a2,c0,a1) T2(w,c0,a0,v0) 
do s_a2 = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_a2,s_a1) .and. & 
IEOR(s_v0,s_a2) == IEOR(s_c0,s_a1) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(v0,a2,c0,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_a1, i_c0) =  &
  V2_(s_a1, s_c0, s_a2)%array(i_a1, i_c0, i_a2)
end do
end do
end do
! Z2 <-- T2(w,c0,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w, i_a0) =  &
  T2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- W9(w,a0,a2,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1))

! W9(w,a0,a2,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W9_(s_a1, s_a2, s_a0, s_w)%array(i_a1, i_a2, i_a0, i_w) = &
    W9_(s_a1, s_a2, s_a0, s_w)%array(i_a1, i_a2, i_a0, i_w) &
  + Z3_(i_a2, i_a1, i_w, i_a0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no3_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no4_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W10, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W10(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xcaaa(sleft, W10, nir, nsym, psym) ! -> Xcaaa (allocate) 
call g_sigma_cooo_ccov_no4_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xcaaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no4_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no4_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W10_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W10_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a2, i_a2, s_a1, i_a1, s_w, i_w, s_a0, i_a0
! W10(w,a0,a2,a1) += (    1.00000000) V2(v0,c0,a2,a1) T2(c0,w,a0,v0) 
do s_c0 = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_a2,s_a1) .and. & 
IEOR(s_v0,s_c0) == IEOR(s_a2,s_a1) .and. &
IEOR(s_c0,s_w) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(v0,c0,a2,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_a1, i_c0) =  &
  V2_(s_a1, s_a2, s_c0)%array(i_a1, i_a2, i_c0)
end do
end do
end do
! Z2 <-- T2(c0,w,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w, i_a0) =  &
  T2_(s_a0, s_w, s_c0)%array(i_a0, i_w, i_c0)
end do
end do
end do

! Z3 <-- W10(w,a0,a2,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1))

! W10(w,a0,a2,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W10_(s_a1, s_a2, s_a0, s_w)%array(i_a1, i_a2, i_a0, i_w) = &
    W10_(s_a1, s_a2, s_a0, s_w)%array(i_a1, i_a2, i_a0, i_w) &
  + Z3_(i_a2, i_a1, i_w, i_a0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no4_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no5_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W11, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W11(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xcaaa(sleft, W11, nir, nsym, psym) ! -> Xcaaa (allocate) 
call g_sigma_cooo_ccov_no5_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xcaaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no5_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no5_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W11_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W11_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a2, i_a2, s_a1, i_a1, s_w, i_w, s_a0, i_a0
! W11(w,a0,a2,a1) += (    1.00000000) V2(v0,c0,a2,a1) T2(w,c0,a0,v0) 
do s_c0 = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_a2,s_a1) .and. & 
IEOR(s_v0,s_c0) == IEOR(s_a2,s_a1) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(v0,c0,a2,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_a1, i_c0) =  &
  V2_(s_a1, s_a2, s_c0)%array(i_a1, i_a2, i_c0)
end do
end do
end do
! Z2 <-- T2(w,c0,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w, i_a0) =  &
  T2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- W11(w,a0,a2,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1))

! W11(w,a0,a2,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W11_(s_a1, s_a2, s_a0, s_w)%array(i_a1, i_a2, i_a0, i_w) = &
    W11_(s_a1, s_a2, s_a0, s_w)%array(i_a1, i_a2, i_a0, i_w) &
  + Z3_(i_a2, i_a1, i_w, i_a0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no5_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no6_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W12, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W12(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xcaaa(sleft, W12, nir, nsym, psym) ! -> Xcaaa (allocate) 
call g_sigma_cooo_ccov_no6_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xcaaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no6_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no6_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W12_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W12_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_i, i_i, s_a1, i_a1, s_w, i_w, s_a0, i_a0
! W12(w,a0,i,a1) += (    1.00000000) V2(v0,c0,i,a1) T2(c0,w,a0,v0) 
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_i,s_a1) .and. & 
IEOR(s_v0,s_c0) == IEOR(s_i,s_a1) .and. &
IEOR(s_c0,s_w) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(v0,c0,i,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a1, i_c0) =  &
  V2_(s_a1, s_i, s_c0)%array(i_a1, i_i, i_c0)
end do
end do
end do
! Z2 <-- T2(c0,w,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w, i_a0) =  &
  T2_(s_a0, s_w, s_c0)%array(i_a0, i_w, i_c0)
end do
end do
end do

! Z3 <-- W12(w,a0,i,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1))

! W12(w,a0,i,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W12_(s_a1, s_i, s_a0, s_w)%array(i_a1, i_i, i_a0, i_w) = &
    W12_(s_a1, s_i, s_a0, s_w)%array(i_a1, i_i, i_a0, i_w) &
  + Z3_(i_i, i_a1, i_w, i_a0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no6_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no7_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W13, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W13(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xcaaa(sleft, W13, nir, nsym, psym) ! -> Xcaaa (allocate) 
call g_sigma_cooo_ccov_no7_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xcaaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no7_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no7_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W13_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W13_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_i, i_i, s_a1, i_a1, s_w, i_w, s_a0, i_a0
! W13(w,a0,i,a1) += (    1.00000000) V2(v0,c0,i,a1) T2(w,c0,a0,v0) 
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_i,s_a1) .and. & 
IEOR(s_v0,s_c0) == IEOR(s_i,s_a1) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(v0,c0,i,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a1, i_c0) =  &
  V2_(s_a1, s_i, s_c0)%array(i_a1, i_i, i_c0)
end do
end do
end do
! Z2 <-- T2(w,c0,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w, i_a0) =  &
  T2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- W13(w,a0,i,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1))

! W13(w,a0,i,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W13_(s_a1, s_i, s_a0, s_w)%array(i_a1, i_i, i_a0, i_w) = &
    W13_(s_a1, s_i, s_a0, s_w)%array(i_a1, i_i, i_a0, i_w) &
  + Z3_(i_i, i_a1, i_w, i_a0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no7_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no8_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W14, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W14(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W14, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no8_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no8_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no8_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W14_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
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

integer :: s_c0, i_c0, s_w, i_w, s_c1, i_c1, s_i, i_i
! W14(w,i) += (    1.00000000) V2(v0,c0,w,c1) T2(c1,c0,i,v0) 
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_c1 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == 0 .and. & 
IEOR(s_v0,s_c0) == IEOR(s_w,s_c1) .and. &
IEOR(s_c1,s_c0) == IEOR(s_i,s_v0)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c1,c0,i,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_c1, i_c0) =  &
  T2_(s_i, s_c0, s_c1)%array(i_i, i_c0, i_c1)
end do
end do
end do
! Z2 <-- V2(v0,c0,w,c1) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z2_(i_c1, i_c0, i_w) =  &
  V2_(s_c1, s_w, s_c0)%array(i_c1, i_w, i_c0)
end do
end do
end do

! Z3 <-- W14(w,i) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! W14(w,i)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W14_(s_i, s_w)%array(i_i, i_w) = &
    W14_(s_i, s_w)%array(i_i, i_w) &
  + Z3_(i_i, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no8_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no9_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W15, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W15(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W15, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no9_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no9_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no9_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W15_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W15_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_w, i_w, s_c1, i_c1, s_i, i_i
! W15(w,i) += (    1.00000000) V2(v0,c0,w,c1) T2(c0,c1,i,v0) 
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_c1 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == 0 .and. & 
IEOR(s_v0,s_c0) == IEOR(s_w,s_c1) .and. &
IEOR(s_c0,s_c1) == IEOR(s_i,s_v0)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1) > 0) then

! Z1 <-- T2(c0,c1,i,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1)))
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_c0, i_c1) =  &
  T2_(s_i, s_c1, s_c0)%array(i_i, i_c1, i_c0)
end do
end do
end do
! Z2 <-- V2(v0,c0,w,c1) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_c1, i_w) =  &
  V2_(s_c1, s_w, s_c0)%array(i_c1, i_w, i_c0)
end do
end do
end do

! Z3 <-- W15(w,i) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! W15(w,i)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W15_(s_i, s_w)%array(i_i, i_w) = &
    W15_(s_i, s_w)%array(i_i, i_w) &
  + Z3_(i_i, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no9_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no10_x0_type1_eri_v &
  (sj, ij, sv0, iv0, V2, W16, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W16(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W16, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no10_x0_type1_eri_v &
  (sj, ij, sv0, iv0, h2_i, Xca, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no10_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no10_x0_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, V2_, W16_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W16_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_c0, i_c0, s_a0, i_a0, s_k, i_k
! W16(c0,k,j,v0) += (    1.00000000) V2(v0,a1,c0,a0) D2(k,j,a1,a0) 
do s_a1 = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0) .and. & 
IEOR(s_v0,s_a1) == IEOR(s_c0,s_a0) .and. &
IEOR(s_k,s_j) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(k,j,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_a1, i_a0) =  &
  D2_(s_a0, s_a1, s_j, s_k)%array(i_a0, i_a1, i_j, i_k)
end do
end do
end do
! Z2 <-- V2(v0,a1,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_c0) =  &
  V2_(s_a0, s_c0, s_a1)%array(i_a0, i_c0, i_a1)
end do
end do
end do

! Z3 <-- W16(c0,k,j,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k))

! W16(c0,k,j,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
W16_(s_k, s_c0)%array(i_k, i_c0) = &
    W16_(s_k, s_c0)%array(i_k, i_c0) &
  + Z3_(i_k, i_c0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no10_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no10_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W16, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W16(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W16, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no10_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no10_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no10_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W16_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
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

integer :: s_w, i_w, s_c0, i_c0, s_i, i_i, s_k, i_k
! S2(w,k,i,j) += (   -2.00000000) T2(w,c0,i,v0) W16(c0,k,j,v0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_c0) == IEOR(s_i,s_v0) .and. &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(w,c0,i,v0) 
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
! Z2 <-- W16(c0,k,j,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_k) =  &
  W16_(s_k, s_c0)%array(i_k, i_c0)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no10_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no11_x0_type1_eri_v &
  (sj, ij, sv0, iv0, V2, W17, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W17(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W17, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no11_x0_type1_eri_v &
  (sj, ij, sv0, iv0, h2_i, Xca, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no11_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no11_x0_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, V2_, W17_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W17_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_c0, i_c0, s_a0, i_a0, s_k, i_k
! W17(c0,k,j,v0) += (    1.00000000) V2(v0,a1,c0,a0) D2(k,j,a1,a0) 
do s_a1 = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0) .and. & 
IEOR(s_v0,s_a1) == IEOR(s_c0,s_a0) .and. &
IEOR(s_k,s_j) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(k,j,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_a1, i_a0) =  &
  D2_(s_a0, s_a1, s_j, s_k)%array(i_a0, i_a1, i_j, i_k)
end do
end do
end do
! Z2 <-- V2(v0,a1,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_c0) =  &
  V2_(s_a0, s_c0, s_a1)%array(i_a0, i_c0, i_a1)
end do
end do
end do

! Z3 <-- W17(c0,k,j,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k))

! W17(c0,k,j,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
W17_(s_k, s_c0)%array(i_k, i_c0) = &
    W17_(s_k, s_c0)%array(i_k, i_c0) &
  + Z3_(i_k, i_c0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no11_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no11_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W17, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W17(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W17, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no11_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no11_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no11_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W17_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
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

integer :: s_c0, i_c0, s_w, i_w, s_i, i_i, s_k, i_k
! S2(w,k,i,j) += (    1.00000000) T2(c0,w,i,v0) W17(c0,k,j,v0) 
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_c0,s_w) == IEOR(s_i,s_v0) .and. &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c0,w,i,v0) 
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
! Z2 <-- W17(c0,k,j,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_k) =  &
  W17_(s_k, s_c0)%array(i_k, i_c0)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no11_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no12_x0_type1_eri_v &
  (sj, ij, sv0, iv0, V2, W18, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W18(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W18, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no12_x0_type1_eri_v &
  (sj, ij, sv0, iv0, h2_i, Xca, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no12_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no12_x0_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, V2_, W18_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W18_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a1, i_a1, s_a0, i_a0, s_k, i_k
! W18(c0,k,j,v0) += (    1.00000000) V2(v0,c0,a1,a0) D2(k,j,a1,a0) 
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0) .and. & 
IEOR(s_v0,s_c0) == IEOR(s_a1,s_a0) .and. &
IEOR(s_k,s_j) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(k,j,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_a1, i_a0) =  &
  D2_(s_a0, s_a1, s_j, s_k)%array(i_a0, i_a1, i_j, i_k)
end do
end do
end do
! Z2 <-- V2(v0,c0,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_c0) =  &
  V2_(s_a0, s_a1, s_c0)%array(i_a0, i_a1, i_c0)
end do
end do
end do

! Z3 <-- W18(c0,k,j,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k))

! W18(c0,k,j,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
W18_(s_k, s_c0)%array(i_k, i_c0) = &
    W18_(s_k, s_c0)%array(i_k, i_c0) &
  + Z3_(i_k, i_c0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no12_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no12_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W18, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W18(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W18, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no12_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no12_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no12_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W18_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W18_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_i, i_i, s_k, i_k
! S2(w,k,i,j) += (    4.00000000) T2(w,c0,i,v0) W18(c0,k,j,v0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_c0) == IEOR(s_i,s_v0) .and. &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(w,c0,i,v0) 
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
! Z2 <-- W18(c0,k,j,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_k) =  &
  W18_(s_k, s_c0)%array(i_k, i_c0)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no12_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no13_x0_type1_eri_v &
  (sj, ij, sv0, iv0, V2, W19, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W19(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W19, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no13_x0_type1_eri_v &
  (sj, ij, sv0, iv0, h2_i, Xca, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no13_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no13_x0_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, V2_, W19_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W19_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a1, i_a1, s_a0, i_a0, s_k, i_k
! W19(c0,k,j,v0) += (    1.00000000) V2(v0,c0,a1,a0) D2(k,j,a1,a0) 
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0) .and. & 
IEOR(s_v0,s_c0) == IEOR(s_a1,s_a0) .and. &
IEOR(s_k,s_j) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(k,j,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_a1, i_a0) =  &
  D2_(s_a0, s_a1, s_j, s_k)%array(i_a0, i_a1, i_j, i_k)
end do
end do
end do
! Z2 <-- V2(v0,c0,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_c0) =  &
  V2_(s_a0, s_a1, s_c0)%array(i_a0, i_a1, i_c0)
end do
end do
end do

! Z3 <-- W19(c0,k,j,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k))

! W19(c0,k,j,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
W19_(s_k, s_c0)%array(i_k, i_c0) = &
    W19_(s_k, s_c0)%array(i_k, i_c0) &
  + Z3_(i_k, i_c0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no13_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no13_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W19, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W19(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W19, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no13_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no13_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no13_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W19_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
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

integer :: s_c0, i_c0, s_w, i_w, s_i, i_i, s_k, i_k
! S2(w,k,i,j) += (   -2.00000000) T2(c0,w,i,v0) W19(c0,k,j,v0) 
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_c0,s_w) == IEOR(s_i,s_v0) .and. &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c0,w,i,v0) 
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
! Z2 <-- W19(c0,k,j,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_k) =  &
  W19_(s_k, s_c0)%array(i_k, i_c0)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no13_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no14_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W20, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W20(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xcaaa(sleft, W20, nir, nsym, psym) ! -> Xcaaa (allocate) 
call g_sigma_cooo_ccov_no14_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xcaaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no14_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no14_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W20_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W20_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_c0, i_c0, s_a1, i_a1, s_w, i_w, s_a0, i_a0
! W20(w,a0,i,a1) += (    1.00000000) V2(v0,i,c0,a1) T2(c0,w,a0,v0) 
do s_i = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_i,s_a1) .and. & 
IEOR(s_v0,s_i) == IEOR(s_c0,s_a1) .and. &
IEOR(s_c0,s_w) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(v0,i,c0,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a1, i_c0) =  &
  V2_(s_a1, s_c0, s_i)%array(i_a1, i_c0, i_i)
end do
end do
end do
! Z2 <-- T2(c0,w,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w, i_a0) =  &
  T2_(s_a0, s_w, s_c0)%array(i_a0, i_w, i_c0)
end do
end do
end do

! Z3 <-- W20(w,a0,i,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1))

! W20(w,a0,i,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W20_(s_a1, s_i, s_a0, s_w)%array(i_a1, i_i, i_a0, i_w) = &
    W20_(s_a1, s_i, s_a0, s_w)%array(i_a1, i_i, i_a0, i_w) &
  + Z3_(i_i, i_a1, i_w, i_a0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no14_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no15_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W21, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W21(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xcaaa(sleft, W21, nir, nsym, psym) ! -> Xcaaa (allocate) 
call g_sigma_cooo_ccov_no15_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xcaaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no15_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no15_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W21_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W21_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_c0, i_c0, s_a1, i_a1, s_w, i_w, s_a0, i_a0
! W21(w,a0,i,a1) += (    1.00000000) V2(v0,i,c0,a1) T2(w,c0,a0,v0) 
do s_i = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_i,s_a1) .and. & 
IEOR(s_v0,s_i) == IEOR(s_c0,s_a1) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(v0,i,c0,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a1, i_c0) =  &
  V2_(s_a1, s_c0, s_i)%array(i_a1, i_c0, i_i)
end do
end do
end do
! Z2 <-- T2(w,c0,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w, i_a0) =  &
  T2_(s_a0, s_c0, s_w)%array(i_a0, i_c0, i_w)
end do
end do
end do

! Z3 <-- W21(w,a0,i,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1))

! W21(w,a0,i,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W21_(s_a1, s_i, s_a0, s_w)%array(i_a1, i_i, i_a0, i_w) = &
    W21_(s_a1, s_i, s_a0, s_w)%array(i_a1, i_i, i_a0, i_w) &
  + Z3_(i_i, i_a1, i_w, i_a0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no15_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no16_x0_type1_eri_v &
  (sv0, iv0, V2, W26, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: V2(*), W26(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W26, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_cooo_ccov_no16_x0_type1_eri_v &
  (sv0, iv0, h2_i, Xcaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no16_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no16_x0_type1_eri_v &
  (s_v0, i_v0, V2_, W26_, D2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W26_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_c0, i_c0, s_a0, i_a0, s_k, i_k, s_i, i_i
! W26(c0,k,i,v0) += (    1.00000000) V2(v0,a1,c0,a0) D2(k,i,a1,a0) 
do s_a1 = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0) .and. & 
IEOR(s_v0,s_a1) == IEOR(s_c0,s_a0) .and. &
IEOR(s_k,s_i) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(k,i,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a1, i_a0) =  &
  D2_(s_a0, s_a1, s_i, s_k)%array(i_a0, i_a1, i_i, i_k)
end do
end do
end do
end do
! Z2 <-- V2(v0,a1,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_c0) =  &
  V2_(s_a0, s_c0, s_a1)%array(i_a0, i_c0, i_a1)
end do
end do
end do

! Z3 <-- W26(c0,k,i,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! W26(c0,k,i,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W26_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) = &
    W26_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) &
  + Z3_(i_k, i_i, i_c0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no16_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no16_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W26, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W26(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W26, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no16_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no16_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no16_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W26_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
integer, intent(in) :: i_v0, s_v0
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
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_k, i_k, s_i, i_i
! S2(w,k,i,j) += (    1.00000000) T2(w,c0,j,v0) W26(c0,k,i,v0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_k = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_c0) == IEOR(s_j,s_v0) .and. &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- W26(c0,k,i,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_c0) =  &
  W26_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0)
end do
end do
end do
! Z2 <-- T2(w,c0,j,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  T2_(s_j, s_c0, s_w)%array(i_j, i_c0, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no16_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no17_x0_type1_eri_v &
  (sv0, iv0, V2, W27, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: V2(*), W27(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W27, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_cooo_ccov_no17_x0_type1_eri_v &
  (sv0, iv0, h2_i, Xcaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no17_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no17_x0_type1_eri_v &
  (s_v0, i_v0, V2_, W27_, D2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W27_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_c0, i_c0, s_a0, i_a0, s_k, i_k, s_i, i_i
! W27(c0,k,i,v0) += (    1.00000000) V2(v0,a1,c0,a0) D2(k,a0,a1,i) 
do s_a1 = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0) .and. & 
IEOR(s_v0,s_a1) == IEOR(s_c0,s_a0) .and. &
IEOR(s_k,s_a0) == IEOR(s_a1,s_i)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(k,a0,a1,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a0, i_a1) =  &
  D2_(s_i, s_a1, s_a0, s_k)%array(i_i, i_a1, i_a0, i_k)
end do
end do
end do
end do
! Z2 <-- V2(v0,a1,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a1, i_c0) =  &
  V2_(s_a0, s_c0, s_a1)%array(i_a0, i_c0, i_a1)
end do
end do
end do

! Z3 <-- W27(c0,k,i,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! W27(c0,k,i,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W27_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) = &
    W27_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) &
  + Z3_(i_k, i_i, i_c0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no17_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no17_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W27, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W27(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sj, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W27, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no17_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no17_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no17_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W27_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
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
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_k, i_k, s_i, i_i
! S2(w,k,i,j) += (    1.00000000) T2(w,c0,v0,j) W27(c0,k,i,v0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_k = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_c0) == IEOR(s_v0,s_j) .and. &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- W27(c0,k,i,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_c0) =  &
  W27_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0)
end do
end do
end do
! Z2 <-- T2(w,c0,v0,j) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  T2_(s_v0, s_c0, s_w)%array(i_v0, i_c0, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no17_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no18_x0_type1_eri_v &
  (sv0, iv0, V2, W28, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: V2(*), W28(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W28, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_cooo_ccov_no18_x0_type1_eri_v &
  (sv0, iv0, h2_i, Xcaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no18_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no18_x0_type1_eri_v &
  (s_v0, i_v0, V2_, W28_, D2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W28_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a1, i_a1, s_a0, i_a0, s_k, i_k, s_i, i_i
! W28(c0,k,i,v0) += (    1.00000000) V2(v0,c0,a1,a0) D2(k,i,a1,a0) 
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0) .and. & 
IEOR(s_v0,s_c0) == IEOR(s_a1,s_a0) .and. &
IEOR(s_k,s_i) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(k,i,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a1, i_a0) =  &
  D2_(s_a0, s_a1, s_i, s_k)%array(i_a0, i_a1, i_i, i_k)
end do
end do
end do
end do
! Z2 <-- V2(v0,c0,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_c0) =  &
  V2_(s_a0, s_a1, s_c0)%array(i_a0, i_a1, i_c0)
end do
end do
end do

! Z3 <-- W28(c0,k,i,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! W28(c0,k,i,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W28_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) = &
    W28_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) &
  + Z3_(i_k, i_i, i_c0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no18_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no18_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W28, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W28(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W28, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no18_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no18_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no18_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W28_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W28_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_k, i_k, s_i, i_i
! S2(w,k,i,j) += (   -2.00000000) T2(w,c0,j,v0) W28(c0,k,i,v0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_k = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_c0) == IEOR(s_j,s_v0) .and. &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- W28(c0,k,i,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_c0) =  &
  W28_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0)
end do
end do
end do
! Z2 <-- T2(w,c0,j,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  T2_(s_j, s_c0, s_w)%array(i_j, i_c0, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no18_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no19_x0_type1_eri_v &
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

call set_symblock_Xcaa(sleft, W29, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_cooo_ccov_no19_x0_type1_eri_v &
  (sv0, iv0, h2_i, Xcaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no19_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no19_x0_type1_eri_v &
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
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a1, i_a1, s_a0, i_a0, s_k, i_k, s_i, i_i
! W29(c0,k,i,v0) += (    1.00000000) V2(v0,c0,a1,a0) D2(k,i,a1,a0) 
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0) .and. & 
IEOR(s_v0,s_c0) == IEOR(s_a1,s_a0) .and. &
IEOR(s_k,s_i) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(k,i,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a1, i_a0) =  &
  D2_(s_a0, s_a1, s_i, s_k)%array(i_a0, i_a1, i_i, i_k)
end do
end do
end do
end do
! Z2 <-- V2(v0,c0,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_c0) =  &
  V2_(s_a0, s_a1, s_c0)%array(i_a0, i_a1, i_c0)
end do
end do
end do

! Z3 <-- W29(c0,k,i,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! W29(c0,k,i,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W29_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) = &
    W29_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) &
  + Z3_(i_k, i_i, i_c0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no19_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no19_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W29, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W29(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sj, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W29, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no19_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no19_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no19_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W29_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
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

integer :: s_w, i_w, s_c0, i_c0, s_k, i_k, s_i, i_i
! S2(w,k,i,j) += (    1.00000000) T2(w,c0,v0,j) W29(c0,k,i,v0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_k = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_c0) == IEOR(s_v0,s_j) .and. &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- W29(c0,k,i,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_c0) =  &
  W29_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0)
end do
end do
end do
! Z2 <-- T2(w,c0,v0,j) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  T2_(s_v0, s_c0, s_w)%array(i_v0, i_c0, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no19_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no20_x0_type1_eri_v &
  (sv0, iv0, V2, W32, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: V2(*), W32(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W32, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_cooo_ccov_no20_x0_type1_eri_v &
  (sv0, iv0, h2_i, Xcaa, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no20_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no20_x0_type1_eri_v &
  (s_v0, i_v0, V2_, W32_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W32_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_i, i_i, s_a0, i_a0, s_k, i_k
! W32(c0,k,i,v0) += (    1.00000000) V2(v0,c0,i,a0) D1(k,a0) 
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0) .and. & 
IEOR(s_v0,s_c0) == IEOR(s_i,s_a0) .and. &
IEOR(s_k,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v0,c0,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z1_(i_c0, i_i, i_a0) =  &
  V2_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0)
end do
end do
end do
! Z2 <-- D1(k,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_k) =  &
  D1_(s_a0, s_k)%array(i_a0, i_k)
end do
end do

! Z3 <-- W32(c0,k,i,v0) 
allocate(Z3_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_i))

! W32(c0,k,i,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W32_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) = &
    W32_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) &
  + Z3_(i_c0, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no20_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no20_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W32, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W32(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W32, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no20_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no20_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no20_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W32_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W32_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_k, i_k, s_i, i_i
! S2(w,k,i,j) += (   -2.00000000) T2(w,c0,j,v0) W32(c0,k,i,v0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_k = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_c0) == IEOR(s_j,s_v0) .and. &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- W32(c0,k,i,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_c0) =  &
  W32_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0)
end do
end do
end do
! Z2 <-- T2(w,c0,j,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  T2_(s_j, s_c0, s_w)%array(i_j, i_c0, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no20_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no21_x0_type1_eri_v &
  (sv0, iv0, V2, W33, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: V2(*), W33(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W33, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_cooo_ccov_no21_x0_type1_eri_v &
  (sv0, iv0, h2_i, Xcaa, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no21_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no21_x0_type1_eri_v &
  (s_v0, i_v0, V2_, W33_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W33_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_i, i_i, s_a0, i_a0, s_k, i_k
! W33(c0,k,i,v0) += (    1.00000000) V2(v0,c0,i,a0) D1(k,a0) 
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0) .and. & 
IEOR(s_v0,s_c0) == IEOR(s_i,s_a0) .and. &
IEOR(s_k,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v0,c0,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z1_(i_c0, i_i, i_a0) =  &
  V2_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0)
end do
end do
end do
! Z2 <-- D1(k,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_k) =  &
  D1_(s_a0, s_k)%array(i_a0, i_k)
end do
end do

! Z3 <-- W33(c0,k,i,v0) 
allocate(Z3_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_i))

! W33(c0,k,i,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W33_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) = &
    W33_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) &
  + Z3_(i_c0, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no21_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no21_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W33, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W33(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sj, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W33, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no21_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no21_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no21_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W33_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W33_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_k, i_k, s_i, i_i
! S2(w,k,i,j) += (    1.00000000) T2(w,c0,v0,j) W33(c0,k,i,v0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_k = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_c0) == IEOR(s_v0,s_j) .and. &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- W33(c0,k,i,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_c0) =  &
  W33_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0)
end do
end do
end do
! Z2 <-- T2(w,c0,v0,j) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  T2_(s_v0, s_c0, s_w)%array(i_v0, i_c0, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no21_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no22_x0_type1_eri_v &
  (sj, ij, sv0, iv0, V2, W34, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W34(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W34, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no22_x0_type1_eri_v &
  (sj, ij, sv0, iv0, h2_i, Xca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no22_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no22_x0_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, V2_, W34_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
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

integer :: s_c0, i_c0, s_a0, i_a0, s_k, i_k
! W34(c0,k,j,v0) += (    1.00000000) V2(v0,c0,j,a0) D1(k,a0) 
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0) .and. & 
IEOR(s_v0,s_c0) == IEOR(s_j,s_a0) .and. &
IEOR(s_k,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D1(k,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_a0) =  &
  D1_(s_a0, s_k)%array(i_a0, i_k)
end do
end do
! Z2 <-- V2(v0,c0,j,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_c0) =  &
  V2_(s_a0, s_j, s_c0)%array(i_a0, i_j, i_c0)
end do
end do

! Z3 <-- W34(c0,k,j,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k))

! W34(c0,k,j,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
W34_(s_k, s_c0)%array(i_k, i_c0) = &
    W34_(s_k, s_c0)%array(i_k, i_c0) &
  + Z3_(i_k, i_c0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no22_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no22_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W34, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W34(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W34, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no22_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no22_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no22_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W34_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W34_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_i, i_i, s_k, i_k
! S2(w,k,i,j) += (    4.00000000) T2(w,c0,i,v0) W34(c0,k,j,v0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_c0) == IEOR(s_i,s_v0) .and. &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(w,c0,i,v0) 
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
! Z2 <-- W34(c0,k,j,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_k) =  &
  W34_(s_k, s_c0)%array(i_k, i_c0)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no22_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no23_x0_type1_eri_v &
  (sj, ij, sv0, iv0, V2, W35, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W35(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W35, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no23_x0_type1_eri_v &
  (sj, ij, sv0, iv0, h2_i, Xca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no23_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no23_x0_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, V2_, W35_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W35_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a0, i_a0, s_k, i_k
! W35(c0,k,j,v0) += (    1.00000000) V2(v0,c0,j,a0) D1(k,a0) 
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0) .and. & 
IEOR(s_v0,s_c0) == IEOR(s_j,s_a0) .and. &
IEOR(s_k,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D1(k,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_a0) =  &
  D1_(s_a0, s_k)%array(i_a0, i_k)
end do
end do
! Z2 <-- V2(v0,c0,j,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_c0) =  &
  V2_(s_a0, s_j, s_c0)%array(i_a0, i_j, i_c0)
end do
end do

! Z3 <-- W35(c0,k,j,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k))

! W35(c0,k,j,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
W35_(s_k, s_c0)%array(i_k, i_c0) = &
    W35_(s_k, s_c0)%array(i_k, i_c0) &
  + Z3_(i_k, i_c0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no23_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no23_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W35, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W35(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W35, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no23_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no23_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no23_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W35_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W35_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_w, i_w, s_i, i_i, s_k, i_k
! S2(w,k,i,j) += (   -2.00000000) T2(c0,w,i,v0) W35(c0,k,j,v0) 
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_c0,s_w) == IEOR(s_i,s_v0) .and. &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c0,w,i,v0) 
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
! Z2 <-- W35(c0,k,j,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_k) =  &
  W35_(s_k, s_c0)%array(i_k, i_c0)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no23_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no24_x0_type1_eri_v &
  (sj, ij, sv0, iv0, V2, W36, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W36(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W36, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no24_x0_type1_eri_v &
  (sj, ij, sv0, iv0, h2_i, Xca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no24_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no24_x0_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, V2_, W36_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W36_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a0, i_a0, s_k, i_k
! W36(c0,k,j,v0) += (    1.00000000) V2(v0,j,c0,a0) D1(k,a0) 
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0) .and. & 
IEOR(s_v0,s_j) == IEOR(s_c0,s_a0) .and. &
IEOR(s_k,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D1(k,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_a0) =  &
  D1_(s_a0, s_k)%array(i_a0, i_k)
end do
end do
! Z2 <-- V2(v0,j,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_c0) =  &
  V2_(s_a0, s_c0, s_j)%array(i_a0, i_c0, i_j)
end do
end do

! Z3 <-- W36(c0,k,j,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k))

! W36(c0,k,j,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
W36_(s_k, s_c0)%array(i_k, i_c0) = &
    W36_(s_k, s_c0)%array(i_k, i_c0) &
  + Z3_(i_k, i_c0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no24_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no24_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W36, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W36(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W36, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no24_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no24_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no24_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W36_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W36_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_i, i_i, s_k, i_k
! S2(w,k,i,j) += (   -2.00000000) T2(w,c0,i,v0) W36(c0,k,j,v0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_c0) == IEOR(s_i,s_v0) .and. &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(w,c0,i,v0) 
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
! Z2 <-- W36(c0,k,j,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_k) =  &
  W36_(s_k, s_c0)%array(i_k, i_c0)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no24_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no25_x0_type1_eri_v &
  (sj, ij, sv0, iv0, V2, W37, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W37(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W37, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no25_x0_type1_eri_v &
  (sj, ij, sv0, iv0, h2_i, Xca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no25_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no25_x0_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, V2_, W37_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W37_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a0, i_a0, s_k, i_k
! W37(c0,k,j,v0) += (    1.00000000) V2(v0,j,c0,a0) D1(k,a0) 
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0) .and. & 
IEOR(s_v0,s_j) == IEOR(s_c0,s_a0) .and. &
IEOR(s_k,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D1(k,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_a0) =  &
  D1_(s_a0, s_k)%array(i_a0, i_k)
end do
end do
! Z2 <-- V2(v0,j,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_c0) =  &
  V2_(s_a0, s_c0, s_j)%array(i_a0, i_c0, i_j)
end do
end do

! Z3 <-- W37(c0,k,j,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k))

! W37(c0,k,j,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
W37_(s_k, s_c0)%array(i_k, i_c0) = &
    W37_(s_k, s_c0)%array(i_k, i_c0) &
  + Z3_(i_k, i_c0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no25_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no25_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W37, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W37(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sv0)

call set_symblock_Xca(sleft, W37, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no25_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xca, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no25_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no25_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W37_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W37_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_w, i_w, s_i, i_i, s_k, i_k
! S2(w,k,i,j) += (    1.00000000) T2(c0,w,i,v0) W37(c0,k,j,v0) 
do s_c0 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_c0,s_w) == IEOR(s_i,s_v0) .and. &
IEOR(s_c0,s_k) == IEOR(s_j,s_v0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c0,w,i,v0) 
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
! Z2 <-- W37(c0,k,j,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_k) =  &
  W37_(s_k, s_c0)%array(i_k, i_c0)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no25_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no26_x0_type1_eri_v &
  (sv0, iv0, V2, W38, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: V2(*), W38(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W38, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_cooo_ccov_no26_x0_type1_eri_v &
  (sv0, iv0, h2_i, Xcaa, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no26_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no26_x0_type1_eri_v &
  (s_v0, i_v0, V2_, W38_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W38_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_c0, i_c0, s_a0, i_a0, s_k, i_k
! W38(c0,k,i,v0) += (    1.00000000) V2(v0,i,c0,a0) D1(k,a0) 
do s_i = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0) .and. & 
IEOR(s_v0,s_i) == IEOR(s_c0,s_a0) .and. &
IEOR(s_k,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v0,i,c0,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_c0, i_a0) =  &
  V2_(s_a0, s_c0, s_i)%array(i_a0, i_c0, i_i)
end do
end do
end do
! Z2 <-- D1(k,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_k) =  &
  D1_(s_a0, s_k)%array(i_a0, i_k)
end do
end do

! Z3 <-- W38(c0,k,i,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_C, s_c0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_C, s_c0))

! W38(c0,k,i,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W38_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) = &
    W38_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) &
  + Z3_(i_i, i_c0, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no26_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no26_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W38, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W38(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sj, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W38, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no26_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no26_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no26_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W38_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W38_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_k, i_k, s_i, i_i
! S2(w,k,i,j) += (   -2.00000000) T2(w,c0,v0,j) W38(c0,k,i,v0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_k = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_c0) == IEOR(s_v0,s_j) .and. &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- W38(c0,k,i,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_c0) =  &
  W38_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0)
end do
end do
end do
! Z2 <-- T2(w,c0,v0,j) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  T2_(s_v0, s_c0, s_w)%array(i_v0, i_c0, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no26_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no27_x0_type1_eri_v &
  (sv0, iv0, V2, W39, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: V2(*), W39(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W39, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_cooo_ccov_no27_x0_type1_eri_v &
  (sv0, iv0, h2_i, Xcaa, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no27_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no27_x0_type1_eri_v &
  (s_v0, i_v0, V2_, W39_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W39_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_c0, i_c0, s_a0, i_a0, s_k, i_k
! W39(c0,k,i,v0) += (    1.00000000) V2(v0,i,c0,a0) D1(k,a0) 
do s_i = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0) .and. & 
IEOR(s_v0,s_i) == IEOR(s_c0,s_a0) .and. &
IEOR(s_k,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v0,i,c0,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_c0, i_a0) =  &
  V2_(s_a0, s_c0, s_i)%array(i_a0, i_c0, i_i)
end do
end do
end do
! Z2 <-- D1(k,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_k) =  &
  D1_(s_a0, s_k)%array(i_a0, i_k)
end do
end do

! Z3 <-- W39(c0,k,i,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_C, s_c0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_C, s_c0))

! W39(c0,k,i,v0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W39_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) = &
    W39_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0) &
  + Z3_(i_i, i_c0, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no27_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no27_x1_type1_eri_v &
  (sj, ij, sv0, iv0, T2, W39, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W39(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W39, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no27_x1_type1_eri_v &
  (sj, ij, sv0, iv0, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_cooo_ccov_no27_x1_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no27_x1_type1_eri_v &
  (s_j, i_j, s_v0, i_v0, T2_, W39_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W39_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_k, i_k, s_i, i_i
! S2(w,k,i,j) += (    1.00000000) T2(w,c0,j,v0) W39(c0,k,i,v0) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_k = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_c0) == IEOR(s_j,s_v0) .and. &
IEOR(s_c0,s_k) == IEOR(s_i,s_v0)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- W39(c0,k,i,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_c0) =  &
  W39_(s_i, s_k, s_c0)%array(i_i, i_k, i_c0)
end do
end do
end do
! Z2 <-- T2(w,c0,j,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_w) =  &
  T2_(s_j, s_c0, s_w)%array(i_j, i_c0, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no27_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no28_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W40, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W40(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W40, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no28_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no28_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no28_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W40_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W40_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_c0, i_c0, s_a0, i_a0, s_w, i_w
! W40(w,a1) += (    1.00000000) V2(v0,a1,c0,a0) T2(w,c0,a0,v0) 
do s_a1 = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_a1) == 0 .and. & 
IEOR(s_v0,s_a1) == IEOR(s_c0,s_a0) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v0,a1,c0,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_c0, i_a0) =  &
  V2_(s_a0, s_c0, s_a1)%array(i_a0, i_c0, i_a1)
end do
end do
end do
! Z2 <-- T2(w,c0,a0,v0) 
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

! Z3 <-- W40(w,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W40(w,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W40_(s_a1, s_w)%array(i_a1, i_w) = &
    W40_(s_a1, s_w)%array(i_a1, i_w) &
  + Z3_(i_a1, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
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

end subroutine g_sigma_cooo_ccov_no28_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no29_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W41, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W41(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W41, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no29_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no29_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no29_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W41_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W41_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_c0, i_c0, s_a0, i_a0, s_w, i_w
! W41(w,a1) += (    1.00000000) V2(v0,a1,c0,a0) T2(c0,w,a0,v0) 
do s_a1 = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_a1) == 0 .and. & 
IEOR(s_v0,s_a1) == IEOR(s_c0,s_a0) .and. &
IEOR(s_c0,s_w) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v0,a1,c0,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_c0, i_a0) =  &
  V2_(s_a0, s_c0, s_a1)%array(i_a0, i_c0, i_a1)
end do
end do
end do
! Z2 <-- T2(c0,w,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0, i_w) =  &
  T2_(s_a0, s_w, s_c0)%array(i_a0, i_w, i_c0)
end do
end do
end do

! Z3 <-- W41(w,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W41(w,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W41_(s_a1, s_w)%array(i_a1, i_w) = &
    W41_(s_a1, s_w)%array(i_a1, i_w) &
  + Z3_(i_a1, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
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

end subroutine g_sigma_cooo_ccov_no29_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no30_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W42, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W42(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W42, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no30_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no30_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no30_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W42_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W42_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a1, i_a1, s_a0, i_a0, s_w, i_w
! W42(w,a1) += (    1.00000000) V2(v0,c0,a1,a0) T2(w,c0,a0,v0) 
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_a1) == 0 .and. & 
IEOR(s_v0,s_c0) == IEOR(s_a1,s_a0) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v0,c0,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_c0, i_a0) =  &
  V2_(s_a0, s_a1, s_c0)%array(i_a0, i_a1, i_c0)
end do
end do
end do
! Z2 <-- T2(w,c0,a0,v0) 
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

! Z3 <-- W42(w,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W42(w,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W42_(s_a1, s_w)%array(i_a1, i_w) = &
    W42_(s_a1, s_w)%array(i_a1, i_w) &
  + Z3_(i_a1, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
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

end subroutine g_sigma_cooo_ccov_no30_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no31_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W43, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W43(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W43, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no31_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no31_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no31_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W43_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W43_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a1, i_a1, s_a0, i_a0, s_w, i_w
! W43(w,a1) += (    1.00000000) V2(v0,c0,a1,a0) T2(c0,w,a0,v0) 
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_a1) == 0 .and. & 
IEOR(s_v0,s_c0) == IEOR(s_a1,s_a0) .and. &
IEOR(s_c0,s_w) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v0,c0,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_c0, i_a0) =  &
  V2_(s_a0, s_a1, s_c0)%array(i_a0, i_a1, i_c0)
end do
end do
end do
! Z2 <-- T2(c0,w,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0, i_w) =  &
  T2_(s_a0, s_w, s_c0)%array(i_a0, i_w, i_c0)
end do
end do
end do

! Z3 <-- W43(w,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W43(w,a1)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W43_(s_a1, s_w)%array(i_a1, i_w) = &
    W43_(s_a1, s_w)%array(i_a1, i_w) &
  + Z3_(i_a1, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
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

end subroutine g_sigma_cooo_ccov_no31_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no32_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W44, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W44(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W44, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no32_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no32_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no32_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W44_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W44_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_i, i_i, s_a0, i_a0, s_w, i_w
! W44(w,i) += (    1.00000000) V2(v0,c0,i,a0) T2(w,c0,a0,v0) 
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_i) == 0 .and. & 
IEOR(s_v0,s_c0) == IEOR(s_i,s_a0) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v0,c0,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_c0, i_a0) =  &
  V2_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0)
end do
end do
end do
! Z2 <-- T2(w,c0,a0,v0) 
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

! Z3 <-- W44(w,i) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! W44(w,i)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W44_(s_i, s_w)%array(i_i, i_w) = &
    W44_(s_i, s_w)%array(i_i, i_w) &
  + Z3_(i_i, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no32_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no33_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W45, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W45(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W45, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no33_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no33_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no33_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W45_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W45_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_i, i_i, s_a0, i_a0, s_w, i_w
! W45(w,i) += (    1.00000000) V2(v0,c0,i,a0) T2(c0,w,a0,v0) 
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_i) == 0 .and. & 
IEOR(s_v0,s_c0) == IEOR(s_i,s_a0) .and. &
IEOR(s_c0,s_w) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v0,c0,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_c0, i_a0) =  &
  V2_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0)
end do
end do
end do
! Z2 <-- T2(c0,w,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0, i_w) =  &
  T2_(s_a0, s_w, s_c0)%array(i_a0, i_w, i_c0)
end do
end do
end do

! Z3 <-- W45(w,i) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! W45(w,i)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W45_(s_i, s_w)%array(i_i, i_w) = &
    W45_(s_i, s_w)%array(i_i, i_w) &
  + Z3_(i_i, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no33_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no34_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W46, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W46(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W46, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no34_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no34_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no34_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W46_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W46_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_c0, i_c0, s_a0, i_a0, s_w, i_w
! W46(w,i) += (    1.00000000) V2(v0,i,c0,a0) T2(c0,w,a0,v0) 
do s_i = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_i) == 0 .and. & 
IEOR(s_v0,s_i) == IEOR(s_c0,s_a0) .and. &
IEOR(s_c0,s_w) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v0,i,c0,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_c0, i_a0) =  &
  V2_(s_a0, s_c0, s_i)%array(i_a0, i_c0, i_i)
end do
end do
end do
! Z2 <-- T2(c0,w,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a0, i_w) =  &
  T2_(s_a0, s_w, s_c0)%array(i_a0, i_w, i_c0)
end do
end do
end do

! Z3 <-- W46(w,i) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! W46(w,i)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W46_(s_i, s_w)%array(i_i, i_w) = &
    W46_(s_i, s_w)%array(i_i, i_w) &
  + Z3_(i_i, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no34_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no35_x0_type1_eri_v &
  (sv0, iv0, T2, V2, W47, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W47(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xca(sleft, W47, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_cooo_ccov_no35_x0_type1_eri_v &
  (sv0, iv0, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no35_x0_type1_eri_v



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
subroutine g_sigma_cooo_ccov_no35_x0_type1_eri_v &
  (s_v0, i_v0, T2_, V2_, W47_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W47_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_c0, i_c0, s_a0, i_a0, s_w, i_w
! W47(w,i) += (    1.00000000) V2(v0,i,c0,a0) T2(w,c0,a0,v0) 
do s_i = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_i) == 0 .and. & 
IEOR(s_v0,s_i) == IEOR(s_c0,s_a0) .and. &
IEOR(s_w,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v0,i,c0,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_c0, i_a0) =  &
  V2_(s_a0, s_c0, s_i)%array(i_a0, i_c0, i_i)
end do
end do
end do
! Z2 <-- T2(w,c0,a0,v0) 
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

! Z3 <-- W47(w,i) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! W47(w,i)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W47_(s_i, s_w)%array(i_i, i_w) = &
    W47_(s_i, s_w)%array(i_i, i_w) &
  + Z3_(i_i, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no35_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no0_x0_type2_eri_v &
  (sj, ij, W6, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W6(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W6, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no0_x0_type2_eri_v &
  (sj, ij, Xca, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no0_x0_type2_eri_v



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
subroutine g_sigma_cooo_ccov_no0_x0_type2_eri_v &
  (s_j, i_j, W6_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W6_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a0, i_a0, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (   -1.00000000) D2(k,j,a0,i) W6(w,a0) 
do s_k = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == IEOR(s_a0,s_i) .and. &
IEOR(s_w,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(k,j,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a0) =  &
  D2_(s_i, s_a0, s_j, s_k)%array(i_i, i_a0, i_j, i_k)
end do
end do
end do
! Z2 <-- W6(w,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  W6_(s_a0, s_w)%array(i_a0, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no0_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no1_x0_type2_eri_v &
  (sj, ij, W7, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W7(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W7, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no1_x0_type2_eri_v &
  (sj, ij, Xca, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no1_x0_type2_eri_v



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
subroutine g_sigma_cooo_ccov_no1_x0_type2_eri_v &
  (s_j, i_j, W7_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W7_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a0, i_a0, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (    2.00000000) D2(k,j,a0,i) W7(w,a0) 
do s_k = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == IEOR(s_a0,s_i) .and. &
IEOR(s_w,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(k,j,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a0) =  &
  D2_(s_i, s_a0, s_j, s_k)%array(i_i, i_a0, i_j, i_k)
end do
end do
end do
! Z2 <-- W7(w,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  W7_(s_a0, s_w)%array(i_a0, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_cooo_ccov_no1_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no2_x0_type2_eri_v &
  (sj, ij, W8, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W8(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xcaaa(sleft, W8, nir, nsym, psym) ! -> Xcaaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no2_x0_type2_eri_v &
  (sj, ij, Xcaaa, av2_i2, d3, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no2_x0_type2_eri_v



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
subroutine g_sigma_cooo_ccov_no2_x0_type2_eri_v &
  (s_j, i_j, W8_, S2_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W8_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a2, i_a2, s_i, i_i, s_a0, i_a0, s_a1, i_a1
integer :: s_w, i_w
! S2(w,k,i,j) += (    1.00000000) D3(k,j,a2,i,a0,a1) W8(w,a0,a2,a1) 
do s_k = 0, nir-1
do s_a2 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(IEOR(s_k,s_j),s_a2) == IEOR(IEOR(s_i,s_a0),s_a1) .and. &
IEOR(s_w,s_a0) == IEOR(s_a2,s_a1)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D3(k,j,a2,i,a0,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a2, i_a0, i_a1) =  &
  D3_(s_a1, s_a0, s_i, s_a2, s_j, s_k)%array(i_a1, i_a0, i_i, i_a2, i_j, i_k)
end do
end do
end do
end do
end do
! Z2 <-- W8(w,a0,a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a0, i_a1, i_w) =  &
  W8_(s_a1, s_a2, s_a0, s_w)%array(i_a1, i_a2, i_a0, i_w)
end do
end do
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no2_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no3_x0_type2_eri_v &
  (sj, ij, W9, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W9(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xcaaa(sleft, W9, nir, nsym, psym) ! -> Xcaaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no3_x0_type2_eri_v &
  (sj, ij, Xcaaa, av2_i2, d3, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no3_x0_type2_eri_v



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
subroutine g_sigma_cooo_ccov_no3_x0_type2_eri_v &
  (s_j, i_j, W9_, S2_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W9_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a2, i_a2, s_a1, i_a1, s_a0, i_a0, s_i, i_i
integer :: s_w, i_w
! S2(w,k,i,j) += (    1.00000000) D3(k,j,a2,a1,a0,i) W9(w,a0,a2,a1) 
do s_k = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(IEOR(s_k,s_j),s_a2) == IEOR(IEOR(s_a1,s_a0),s_i) .and. &
IEOR(s_w,s_a0) == IEOR(s_a2,s_a1)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D3(k,j,a2,a1,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a2, i_a1, i_a0) =  &
  D3_(s_i, s_a0, s_a1, s_a2, s_j, s_k)%array(i_i, i_a0, i_a1, i_a2, i_j, i_k)
end do
end do
end do
end do
end do
! Z2 <-- W9(w,a0,a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1, i_a0, i_w) =  &
  W9_(s_a1, s_a2, s_a0, s_w)%array(i_a1, i_a2, i_a0, i_w)
end do
end do
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no3_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no4_x0_type2_eri_v &
  (sj, ij, W10, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W10(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xcaaa(sleft, W10, nir, nsym, psym) ! -> Xcaaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no4_x0_type2_eri_v &
  (sj, ij, Xcaaa, av2_i2, d3, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no4_x0_type2_eri_v



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
subroutine g_sigma_cooo_ccov_no4_x0_type2_eri_v &
  (s_j, i_j, W10_, S2_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W10_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a2, i_a2, s_a1, i_a1, s_a0, i_a0, s_i, i_i
integer :: s_w, i_w
! S2(w,k,i,j) += (    1.00000000) D3(k,j,a2,a1,a0,i) W10(w,a0,a2,a1) 
do s_k = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(IEOR(s_k,s_j),s_a2) == IEOR(IEOR(s_a1,s_a0),s_i) .and. &
IEOR(s_w,s_a0) == IEOR(s_a2,s_a1)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D3(k,j,a2,a1,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a2, i_a1, i_a0) =  &
  D3_(s_i, s_a0, s_a1, s_a2, s_j, s_k)%array(i_i, i_a0, i_a1, i_a2, i_j, i_k)
end do
end do
end do
end do
end do
! Z2 <-- W10(w,a0,a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1, i_a0, i_w) =  &
  W10_(s_a1, s_a2, s_a0, s_w)%array(i_a1, i_a2, i_a0, i_w)
end do
end do
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no4_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no5_x0_type2_eri_v &
  (sj, ij, W11, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W11(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xcaaa(sleft, W11, nir, nsym, psym) ! -> Xcaaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no5_x0_type2_eri_v &
  (sj, ij, Xcaaa, av2_i2, d3, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no5_x0_type2_eri_v



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
subroutine g_sigma_cooo_ccov_no5_x0_type2_eri_v &
  (s_j, i_j, W11_, S2_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W11_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a2, i_a2, s_a1, i_a1, s_a0, i_a0, s_i, i_i
integer :: s_w, i_w
! S2(w,k,i,j) += (   -2.00000000) D3(k,j,a2,a1,a0,i) W11(w,a0,a2,a1) 
do s_k = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(IEOR(s_k,s_j),s_a2) == IEOR(IEOR(s_a1,s_a0),s_i) .and. &
IEOR(s_w,s_a0) == IEOR(s_a2,s_a1)) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D3(k,j,a2,a1,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a2, i_a1, i_a0) =  &
  D3_(s_i, s_a0, s_a1, s_a2, s_j, s_k)%array(i_i, i_a0, i_a1, i_a2, i_j, i_k)
end do
end do
end do
end do
end do
! Z2 <-- W11(w,a0,a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1, i_a0, i_w) =  &
  W11_(s_a1, s_a2, s_a0, s_w)%array(i_a1, i_a2, i_a0, i_w)
end do
end do
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no5_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no6_x0_type2_eri_v &
  (sj, ij, W12, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W12(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xcaaa(sleft, W12, nir, nsym, psym) ! -> Xcaaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no6_x0_type2_eri_v &
  (sj, ij, Xcaaa, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no6_x0_type2_eri_v



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
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_cooo_ccov_no6_x0_type2_eri_v &
  (s_j, i_j, W12_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W12_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a0, i_a0, s_a1, i_a1, s_w, i_w, s_i, i_i
! S2(w,k,i,j) += (    1.00000000) D2(k,j,a0,a1) W12(w,a0,i,a1) 
do s_k = 0, nir-1
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == IEOR(s_a0,s_a1) .and. &
IEOR(s_w,s_a0) == IEOR(s_i,s_a1)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- W12(w,a0,i,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i, i_a0, i_a1) =  &
  W12_(s_a1, s_i, s_a0, s_w)%array(i_a1, i_i, i_a0, i_w)
end do
end do
end do
end do
! Z2 <-- D2(k,j,a0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a1, i_k) =  &
  D2_(s_a1, s_a0, s_j, s_k)%array(i_a1, i_a0, i_j, i_k)
end do
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no6_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no7_x0_type2_eri_v &
  (sj, ij, W13, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W13(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xcaaa(sleft, W13, nir, nsym, psym) ! -> Xcaaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no7_x0_type2_eri_v &
  (sj, ij, Xcaaa, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no7_x0_type2_eri_v



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
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_cooo_ccov_no7_x0_type2_eri_v &
  (s_j, i_j, W13_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W13_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a0, i_a0, s_a1, i_a1, s_w, i_w, s_i, i_i
! S2(w,k,i,j) += (   -2.00000000) D2(k,j,a0,a1) W13(w,a0,i,a1) 
do s_k = 0, nir-1
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == IEOR(s_a0,s_a1) .and. &
IEOR(s_w,s_a0) == IEOR(s_i,s_a1)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- W13(w,a0,i,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i, i_a0, i_a1) =  &
  W13_(s_a1, s_i, s_a0, s_w)%array(i_a1, i_i, i_a0, i_w)
end do
end do
end do
end do
! Z2 <-- D2(k,j,a0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a1, i_k) =  &
  D2_(s_a1, s_a0, s_j, s_k)%array(i_a1, i_a0, i_j, i_k)
end do
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no7_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no8_x0_type2_eri_v &
  (sj, ij, W14, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W14(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W14, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no8_x0_type2_eri_v &
  (sj, ij, Xca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no8_x0_type2_eri_v



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
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_cooo_ccov_no8_x0_type2_eri_v &
  (s_j, i_j, W14_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
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
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_w, i_w, s_i, i_i
! S2(w,k,i,j) += (   -4.00000000) D1(k,j) W14(w,i) 
do s_k = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == 0 .and. &
IEOR(s_w,s_i) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0) then

! Z1 <-- W14(w,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i) =  &
  W14_(s_i, s_w)%array(i_i, i_w)
end do
end do
! Z2 <-- D1(k,j) 
allocate(Z2_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_k) =  &
  D1_(s_j, s_k)%array(i_j, i_k)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     1,&
                     - 4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no8_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no9_x0_type2_eri_v &
  (sj, ij, W15, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W15(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W15, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no9_x0_type2_eri_v &
  (sj, ij, Xca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no9_x0_type2_eri_v



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
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_cooo_ccov_no9_x0_type2_eri_v &
  (s_j, i_j, W15_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W15_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_w, i_w, s_i, i_i
! S2(w,k,i,j) += (    2.00000000) D1(k,j) W15(w,i) 
do s_k = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == 0 .and. &
IEOR(s_w,s_i) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0) then

! Z1 <-- W15(w,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i) =  &
  W15_(s_i, s_w)%array(i_i, i_w)
end do
end do
! Z2 <-- D1(k,j) 
allocate(Z2_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_k) =  &
  D1_(s_j, s_k)%array(i_j, i_k)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     1,&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no9_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no10_x0_type2_eri_v &
  (sj, ij, W20, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W20(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xcaaa(sleft, W20, nir, nsym, psym) ! -> Xcaaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no10_x0_type2_eri_v &
  (sj, ij, Xcaaa, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no10_x0_type2_eri_v



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
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_cooo_ccov_no10_x0_type2_eri_v &
  (s_j, i_j, W20_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W20_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a0, i_a0, s_a1, i_a1, s_w, i_w, s_i, i_i
! S2(w,k,i,j) += (   -2.00000000) D2(k,j,a0,a1) W20(w,a0,i,a1) 
do s_k = 0, nir-1
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == IEOR(s_a0,s_a1) .and. &
IEOR(s_w,s_a0) == IEOR(s_i,s_a1)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- W20(w,a0,i,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i, i_a0, i_a1) =  &
  W20_(s_a1, s_i, s_a0, s_w)%array(i_a1, i_i, i_a0, i_w)
end do
end do
end do
end do
! Z2 <-- D2(k,j,a0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a1, i_k) =  &
  D2_(s_a1, s_a0, s_j, s_k)%array(i_a1, i_a0, i_j, i_k)
end do
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no10_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no11_x0_type2_eri_v &
  (sj, ij, W21, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W21(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xcaaa(sleft, W21, nir, nsym, psym) ! -> Xcaaa (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no11_x0_type2_eri_v &
  (sj, ij, Xcaaa, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaaa)

end subroutine g_if_sigma_cooo_ccov_no11_x0_type2_eri_v



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
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_cooo_ccov_no11_x0_type2_eri_v &
  (s_j, i_j, W21_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W21_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a0, i_a0, s_a1, i_a1, s_w, i_w, s_i, i_i
! S2(w,k,i,j) += (    1.00000000) D2(k,j,a0,a1) W21(w,a0,i,a1) 
do s_k = 0, nir-1
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == IEOR(s_a0,s_a1) .and. &
IEOR(s_w,s_a0) == IEOR(s_i,s_a1)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- W21(w,a0,i,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i, i_a0, i_a1) =  &
  W21_(s_a1, s_i, s_a0, s_w)%array(i_a1, i_i, i_a0, i_w)
end do
end do
end do
end do
! Z2 <-- D2(k,j,a0,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a1, i_k) =  &
  D2_(s_a1, s_a0, s_j, s_k)%array(i_a1, i_a0, i_j, i_k)
end do
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no11_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no12_x0_type2_eri_v &
  (sj, ij, W40, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W40(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W40, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no12_x0_type2_eri_v &
  (sj, ij, Xca, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no12_x0_type2_eri_v



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
subroutine g_sigma_cooo_ccov_no12_x0_type2_eri_v &
  (s_j, i_j, W40_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W40_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a1, i_a1, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (    1.00000000) D2(k,j,a1,i) W40(w,a1) 
do s_k = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == IEOR(s_a1,s_i) .and. &
IEOR(s_w,s_a1) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(k,j,a1,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a1) =  &
  D2_(s_i, s_a1, s_j, s_k)%array(i_i, i_a1, i_j, i_k)
end do
end do
end do
! Z2 <-- W40(w,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w) =  &
  W40_(s_a1, s_w)%array(i_a1, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no12_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no13_x0_type2_eri_v &
  (sj, ij, W41, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W41(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W41, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no13_x0_type2_eri_v &
  (sj, ij, Xca, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no13_x0_type2_eri_v



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
subroutine g_sigma_cooo_ccov_no13_x0_type2_eri_v &
  (s_j, i_j, W41_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W41_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a1, i_a1, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (   -2.00000000) D2(k,j,a1,i) W41(w,a1) 
do s_k = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == IEOR(s_a1,s_i) .and. &
IEOR(s_w,s_a1) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(k,j,a1,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a1) =  &
  D2_(s_i, s_a1, s_j, s_k)%array(i_i, i_a1, i_j, i_k)
end do
end do
end do
! Z2 <-- W41(w,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w) =  &
  W41_(s_a1, s_w)%array(i_a1, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no13_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no14_x0_type2_eri_v &
  (sj, ij, W42, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W42(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W42, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no14_x0_type2_eri_v &
  (sj, ij, Xca, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no14_x0_type2_eri_v



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
subroutine g_sigma_cooo_ccov_no14_x0_type2_eri_v &
  (s_j, i_j, W42_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W42_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a1, i_a1, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (   -2.00000000) D2(k,j,a1,i) W42(w,a1) 
do s_k = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == IEOR(s_a1,s_i) .and. &
IEOR(s_w,s_a1) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(k,j,a1,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a1) =  &
  D2_(s_i, s_a1, s_j, s_k)%array(i_i, i_a1, i_j, i_k)
end do
end do
end do
! Z2 <-- W42(w,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w) =  &
  W42_(s_a1, s_w)%array(i_a1, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no14_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no15_x0_type2_eri_v &
  (sj, ij, W43, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W43(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W43, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no15_x0_type2_eri_v &
  (sj, ij, Xca, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no15_x0_type2_eri_v



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
subroutine g_sigma_cooo_ccov_no15_x0_type2_eri_v &
  (s_j, i_j, W43_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W43_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_a1, i_a1, s_i, i_i, s_w, i_w
! S2(w,k,i,j) += (    1.00000000) D2(k,j,a1,i) W43(w,a1) 
do s_k = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == IEOR(s_a1,s_i) .and. &
IEOR(s_w,s_a1) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(k,j,a1,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_a1) =  &
  D2_(s_i, s_a1, s_j, s_k)%array(i_i, i_a1, i_j, i_k)
end do
end do
end do
! Z2 <-- W43(w,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w) =  &
  W43_(s_a1, s_w)%array(i_a1, i_w)
end do
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_k, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no15_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no16_x0_type2_eri_v &
  (sj, ij, W44, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W44(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W44, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no16_x0_type2_eri_v &
  (sj, ij, Xca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no16_x0_type2_eri_v



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
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_cooo_ccov_no16_x0_type2_eri_v &
  (s_j, i_j, W44_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W44_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_w, i_w, s_i, i_i
! S2(w,k,i,j) += (    4.00000000) D1(k,j) W44(w,i) 
do s_k = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == 0 .and. &
IEOR(s_w,s_i) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0) then

! Z1 <-- W44(w,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i) =  &
  W44_(s_i, s_w)%array(i_i, i_w)
end do
end do
! Z2 <-- D1(k,j) 
allocate(Z2_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_k) =  &
  D1_(s_j, s_k)%array(i_j, i_k)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     1,&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no16_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no17_x0_type2_eri_v &
  (sj, ij, W45, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W45(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W45, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no17_x0_type2_eri_v &
  (sj, ij, Xca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no17_x0_type2_eri_v



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
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_cooo_ccov_no17_x0_type2_eri_v &
  (s_j, i_j, W45_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W45_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_w, i_w, s_i, i_i
! S2(w,k,i,j) += (   -2.00000000) D1(k,j) W45(w,i) 
do s_k = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == 0 .and. &
IEOR(s_w,s_i) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0) then

! Z1 <-- W45(w,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i) =  &
  W45_(s_i, s_w)%array(i_i, i_w)
end do
end do
! Z2 <-- D1(k,j) 
allocate(Z2_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_k) =  &
  D1_(s_j, s_k)%array(i_j, i_k)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     1,&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no17_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no18_x0_type2_eri_v &
  (sj, ij, W46, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W46(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W46, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no18_x0_type2_eri_v &
  (sj, ij, Xca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no18_x0_type2_eri_v



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
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_cooo_ccov_no18_x0_type2_eri_v &
  (s_j, i_j, W46_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W46_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_w, i_w, s_i, i_i
! S2(w,k,i,j) += (    4.00000000) D1(k,j) W46(w,i) 
do s_k = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == 0 .and. &
IEOR(s_w,s_i) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0) then

! Z1 <-- W46(w,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i) =  &
  W46_(s_i, s_w)%array(i_i, i_w)
end do
end do
! Z2 <-- D1(k,j) 
allocate(Z2_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_k) =  &
  D1_(s_j, s_k)%array(i_j, i_k)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     1,&
                     4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no18_x0_type2_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_cooo_ccov_no19_x0_type2_eri_v &
  (sj, ij, W47, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: W47(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xca(sleft, W47, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_cooo_ccov_no19_x0_type2_eri_v &
  (sj, ij, Xca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_cooo_ccov_no19_x0_type2_eri_v



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
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_cooo_ccov_no19_x0_type2_eri_v &
  (s_j, i_j, W47_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W47_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_w, i_w, s_i, i_i
! S2(w,k,i,j) += (   -2.00000000) D1(k,j) W47(w,i) 
do s_k = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_k) == IEOR(s_i,s_j) .and. & 
IEOR(s_k,s_j) == 0 .and. &
IEOR(s_w,s_i) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_k) > 0) then

! Z1 <-- W47(w,i) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i) =  &
  W47_(s_i, s_w)%array(i_i, i_w)
end do
end do
! Z2 <-- D1(k,j) 
allocate(Z2_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_k) =  &
  D1_(s_j, s_k)%array(i_j, i_k)
end do

! Z3 <-- S2(w,k,i,j) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_k),&
                     1,&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i))

! S2(w,k,i,j)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) = &
    S2_(s_i, s_k, s_w)%array(i_i, i_k, i_w) &
  + Z3_(i_w, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_k) * &
                1 * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_cooo_ccov_no19_x0_type2_eri_v

