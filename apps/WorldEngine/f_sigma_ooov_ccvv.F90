#include <sci/icmr/fsrc/f_mr.fh>


!    o__ __o__/_                            o                         
!   <|    v                                <|>                        
!   < >                                    < >                        
!    |         o__  __o   \o__ __o__ __o    |        o__ __o         
!    o__/_    /v      |>   |     |     |>   o__/_   /v     v\        
!    |       />      //   / \   / \   / \   |      />       <\    
!   <o>      \o    o/     \o/   \o/   \o/   |      \         /   
!    |        v\  /v __o   |     |     |    o       o       o        
!   / \        <\/> __/>  / \   / \   / \   <\__    <\__ __/>  

!                                    Generated date : Sun Apr 20 10:26:09 2014



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ccvv_no0_x0_type1_eri_o &
  (sa, ia, sa0, ia0, T2, V2, W0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), W0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call g_sigma_ooov_ccvv_no0_x0_type1_eri_o &
  (sa, ia, sa0, ia0, av2_i, h2_i, W0, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_ooov_ccvv_no0_x0_type1_eri_o



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
subroutine g_sigma_ooov_ccvv_no0_x0_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, W0_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W0_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_v0, i_v0, s_c0, i_c0
! W0(a0,a) += (    1.00000000) V2(a0,c1,v0,c0) T2(c0,c1,v0,a) 
do s_c1 = 0, nir-1
do s_v0 = 0, nir-1
do s_c0 = 0, nir-1
if( &
IEOR(s_a0,s_a) == 0 .and. & 
IEOR(s_a0,s_c1) == IEOR(s_v0,s_c0) .and. &
IEOR(s_c0,s_c1) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(a0,c1,v0,c0) 
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
! Z2 <-- T2(c0,c1,v0,a) 
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

! W0(a0,a)  <-- Z3
W0_ = &
    W0_ &
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

end subroutine g_sigma_ooov_ccvv_no0_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ccvv_no0_x1_type1_eri_o &
  (sa, ia, sa0, ia0, W0, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: W0, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa0,sa)

call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ccvv_no0_x1_type1_eri_o &
  (sa, ia, sa0, ia0, W0, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ccvv_no0_x1_type1_eri_o



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
subroutine g_sigma_ooov_ccvv_no0_x1_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, W0_, S2_, D2_, nir, nsym, psym, flops)

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
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W0_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_i, i_i, s_j, i_j
! S2(i,j,k,a) += (   -4.00000000) D2(k,i,a0,j) W0(a0,a) 
do s_k = 0, nir-1
do s_i = 0, nir-1
do s_j = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_k,s_a) .and. & 
IEOR(s_k,s_i) == IEOR(s_a0,s_j) .and. &
IEOR(s_a0,s_a) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j) > 0) then

! Z1 <-- D2(k,i,a0,j) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j)))
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_j) =  &
  D2_(s_j, s_a0, s_i, s_k)%array(i_j, i_a0, i_i, i_k)
end do
end do
end do
! Z2 <-- W0(a0,a) 
Z2_ =  &
  W0_

! Z3 <-- S2(i,j,k,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j),&
                     1,&
                     1,&
                     - 4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j))

! S2(i,j,k,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_k, s_j, s_i)%array(i_k, i_j, i_i) = &
    S2_(s_k, s_j, s_i)%array(i_k, i_j, i_i) &
  + Z3_(i_k, i_i, i_j)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j) * &
                1 * &
                1 * 2.0d+00

deallocate(Z1_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ccvv_no0_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ccvv_no1_x0_type1_eri_o &
  (sa, ia, sa0, ia0, T2, V2, W1, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), W1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call g_sigma_ooov_ccvv_no1_x0_type1_eri_o &
  (sa, ia, sa0, ia0, av2_i, h2_i, W1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_ooov_ccvv_no1_x0_type1_eri_o



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
subroutine g_sigma_ooov_ccvv_no1_x0_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, W1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W1_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_v0, i_v0, s_c0, i_c0
! W1(a0,a) += (    1.00000000) V2(a0,c1,v0,c0) T2(c1,c0,v0,a) 
do s_c1 = 0, nir-1
do s_v0 = 0, nir-1
do s_c0 = 0, nir-1
if( &
IEOR(s_a0,s_a) == 0 .and. & 
IEOR(s_a0,s_c1) == IEOR(s_v0,s_c0) .and. &
IEOR(s_c1,s_c0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_V, s_v0)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(a0,c1,v0,c0) 
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
! Z2 <-- T2(c1,c0,v0,a) 
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

! W1(a0,a)  <-- Z3
W1_ = &
    W1_ &
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

end subroutine g_sigma_ooov_ccvv_no1_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ccvv_no1_x1_type1_eri_o &
  (sa, ia, sa0, ia0, W1, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: W1, S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa0,sa)

call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ccvv_no1_x1_type1_eri_o &
  (sa, ia, sa0, ia0, W1, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ccvv_no1_x1_type1_eri_o



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
subroutine g_sigma_ooov_ccvv_no1_x1_type1_eri_o &
  (s_a, i_a, s_a0, i_a0, W1_, S2_, D2_, nir, nsym, psym, flops)

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
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: W1_

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_i, i_i, s_j, i_j
! S2(i,j,k,a) += (    2.00000000) D2(k,i,a0,j) W1(a0,a) 
do s_k = 0, nir-1
do s_i = 0, nir-1
do s_j = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_k,s_a) .and. & 
IEOR(s_k,s_i) == IEOR(s_a0,s_j) .and. &
IEOR(s_a0,s_a) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j) > 0) then

! Z1 <-- D2(k,i,a0,j) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j)))
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_j) =  &
  D2_(s_j, s_a0, s_i, s_k)%array(i_j, i_a0, i_i, i_k)
end do
end do
end do
! Z2 <-- W1(a0,a) 
Z2_ =  &
  W1_

! Z3 <-- S2(i,j,k,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j),&
                     1,&
                     1,&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j))

! S2(i,j,k,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_k, s_j, s_i)%array(i_k, i_j, i_i) = &
    S2_(s_k, s_j, s_i)%array(i_k, i_j, i_i) &
  + Z3_(i_k, i_i, i_j)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j) * &
                1 * &
                1 * 2.0d+00

deallocate(Z1_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_ccvv_no1_x1_type1_eri_o

