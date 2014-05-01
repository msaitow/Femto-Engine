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
subroutine g_if_sigma_ooov_g_no0_x0_type0_noeri &
  (sa, ia, T0, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T0
real(kind=8), intent(inout) :: S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_g_no0_x0_type0_noeri &
  (sa, ia, T0, av2_i2, d2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)

end subroutine g_if_sigma_ooov_g_no0_x0_type0_noeri



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
subroutine g_sigma_ooov_g_no0_x0_type0_noeri &
  (s_a, i_a, T0, S2_, D2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: T0
! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_i, i_i, s_a0, i_a0, s_j, i_j
! S2(i,j,k,a) += (    1.00000000) T0 D2(k,i,a0,j) Fc1(a0,a) 
do s_k = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_j = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_k,s_a) .and. & 
IEOR(s_k,s_i) == IEOR(s_a0,s_j) .and. &
IEOR(s_a0,s_a) == 0) then

if(psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(k,i,a0,j) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_i, i_j, i_a0) =  &
  D2_(s_j, s_a0, s_i, s_k)%array(i_j, i_a0, i_i, i_k)
end do
end do
end do
end do
! Z2 <-- Fc1(a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0) =  &
  Fc1_(s_a, s_a0)%array(i_a, i_a0)
end do

! Z3 <-- S2(i,j,k,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j),&
                     1,&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00*T0, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
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
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_g_no0_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_g_no0_x0_type1_eri_v &
  (sa, ia, T0, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T0
real(kind=8), intent(inout) :: V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_g_no0_x0_type1_eri_v &
  (sa, ia, T0, h2_i, av2_i2, d3, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)

end subroutine g_if_sigma_ooov_g_no0_x0_type1_eri_v



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Information for debugging 
! >> Score :: 111111
! The best shot!
! RDM is rotated :: D3(k,i,a2,j,a1,a0)  >> D3(a1,a0,a2,j,k,i) 
! rowInd : @[i, "active"] @[k, "active"] @[j, "active"] 
! summedInd : @[a2, "active"] @[a0, "active"] @[a1, "active"] 
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_ooov_g_no0_x0_type1_eri_v &
  (s_a, i_a, T0, V2_, S2_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: T0
! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a1, i_a1, s_a0, i_a0, s_k, i_k, s_i, i_i
integer :: s_j, i_j
! S2(i,j,k,a) += (    1.00000000) T0 V2(a,a2,a1,a0) D3(k,i,a2,j,a1,a0) 
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_k = 0, nir-1
do s_i = 0, nir-1
do s_j = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_k,s_a) .and. & 
IEOR(s_a,s_a2) == IEOR(s_a1,s_a0) .and. &
IEOR(IEOR(s_k,s_i),s_a2) == IEOR(IEOR(s_j,s_a1),s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z2 <-- V2(a,a2,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a0, i_a1) =  &
  V2_(s_a0, s_a1, s_a2)%array(i_a0, i_a1, i_a2)
end do
end do
end do

! Z3 <-- S2(i,j,k,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_j),&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00*T0, &
                     D3_(s_i, s_k, s_j, s_a2, s_a0, s_a1)%array,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_j))

! S2(i,j,k,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_k, s_j, s_i)%array(i_k, i_j, i_i) = &
    S2_(s_k, s_j, s_i)%array(i_k, i_j, i_i) &
  + Z3_(i_i, i_k, i_j)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_j) * &
                1 * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ooov_g_no0_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_g_no1_x0_type1_eri_v &
  (sa, ia, T0, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T0
real(kind=8), intent(inout) :: V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_g_no1_x0_type1_eri_v &
  (sa, ia, T0, h2_i, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)

end subroutine g_if_sigma_ooov_g_no1_x0_type1_eri_v



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
subroutine g_sigma_ooov_g_no1_x0_type1_eri_v &
  (s_a, i_a, T0, V2_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: T0
! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_k, i_k, s_a0, i_a0, s_j, i_j, s_i, i_i
! S2(i,j,k,a) += (    1.00000000) T0 V2(a,a1,k,a0) D2(j,a1,i,a0) 
do s_a1 = 0, nir-1
do s_k = 0, nir-1
do s_a0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_k,s_a) .and. & 
IEOR(s_a,s_a1) == IEOR(s_k,s_a0) .and. &
IEOR(s_j,s_a1) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(a,a1,k,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z1_(i_k, i_a1, i_a0) =  &
  V2_(s_a0, s_k, s_a1)%array(i_a0, i_k, i_a1)
end do
end do
end do
! Z2 <-- D2(j,a1,i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_j, i_i) =  &
  D2_(s_a0, s_i, s_a1, s_j)%array(i_a0, i_i, i_a1, i_j)
end do
end do
end do
end do

! Z3 <-- S2(i,j,k,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00*T0, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_k),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_k))

! S2(i,j,k,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
S2_(s_k, s_j, s_i)%array(i_k, i_j, i_i) = &
    S2_(s_k, s_j, s_i)%array(i_k, i_j, i_i) &
  + Z3_(i_k, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_ooov_g_no1_x0_type1_eri_v

