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
subroutine g_if_sigma_oovv_ooov_no0_x0_type0_noeri &
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

call set_symblock_Xaaa(sleft, W0, nir, nsym, psym) ! -> Xaaa (allocate) 
call g_sigma_oovv_ooov_no0_x0_type0_noeri &
  (sb, ib, av2_i, Xaaa, d3, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x0_type0_noeri



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Information for debugging 
! >> Score :: 111111
! The best shot!
! RDM is rotated :: D3(j,a1,i,a3,a2,a0)  >> D3(a0,a2,a1,j,a3,i) 
! rowInd : @[i, "active"] @[a3, "active"] @[j, "active"] 
! summedInd : @[a1, "active"] @[a2, "active"] @[a0, "active"] 
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_oovv_ooov_no0_x0_type0_noeri &
  (s_b, i_b, T2_, W0_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_a1, i_a1, s_a2, i_a2, s_j, i_j, s_i, i_i
integer :: s_a3, i_a3
! W0(j,i,a3,b) += (    1.00000000) T2(a0,a1,a2,b) D3(j,a1,i,a3,a2,a0) 
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a2 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a3 = 0, nir-1
if( &
IEOR(s_j,s_i) == IEOR(s_a3,s_b) .and. & 
IEOR(s_a0,s_a1) == IEOR(s_a2,s_b) .and. &
IEOR(IEOR(s_j,s_a1),s_i) == IEOR(IEOR(s_a3,s_a2),s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z2 <-- T2(a0,a1,a2,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a2, i_a0) =  &
  T2_(s_a2, s_a1, s_a0)%array(i_a2, i_a1, i_a0)
end do
end do
end do

! Z3 <-- W0(j,i,a3,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_j),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     D3_(s_i, s_a3, s_j, s_a1, s_a2, s_a0)%array,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_j))

! W0(j,i,a3,b)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
W0_(s_a3, s_i, s_j)%array(i_a3, i_i, i_j) = &
    W0_(s_a3, s_i, s_j)%array(i_a3, i_i, i_j) &
  + Z3_(i_i, i_a3, i_j)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_j) * &
                1 * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

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

end subroutine g_sigma_oovv_ooov_no0_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x1_type0_noeri &
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

call set_symblock_Xaaa(sleft, W0, nir, nsym, psym) ! -> Xaaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no0_x1_type0_noeri &
  (sb, ib, Xaaa, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x1_type0_noeri



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
subroutine g_sigma_oovv_ooov_no0_x1_type0_noeri &
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
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_a, i_a, s_j, i_j, s_i, i_i
! S2(i,j,a,b) += (    1.00000000) Fc1(a3,a) W0(j,i,a3,b) 
do s_a3 = 0, nir-1
do s_a = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_a3,s_a) == 0 .and. &
IEOR(s_j,s_i) == IEOR(s_a3,s_b)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3) > 0) then

! Z1 <-- Fc1(a3,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a3) =  &
  Fc1_(s_a, s_a3)%array(i_a, i_a3)
end do
end do
! Z2 <-- W0(j,i,a3,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_j, i_i) =  &
  W0_(s_a3, s_i, s_j)%array(i_a3, i_i, i_j)
end do
end do
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a3),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,j,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) = &
    S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) &
  + Z3_(i_a, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a3) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no0_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x0_type0_noeri &
  (sa2, ia2, T2, W1, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa2, ia2
real(kind=8), intent(inout) :: T2(*), W1(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa2, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xaaav(sleft, W1, nir, nsym, psym) ! -> Xaaav (allocate) 
call g_sigma_oovv_ooov_no1_x0_type0_noeri &
  (sa2, ia2, av2_i, Xaaav, d3, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaaav)

end subroutine g_if_sigma_oovv_ooov_no1_x0_type0_noeri



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
subroutine g_sigma_oovv_ooov_no1_x0_type0_noeri &
  (s_a2, i_a2, T2_, W1_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a2, s_a2
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W1_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a0, i_a0, s_a, i_a, s_j, i_j, s_a3, i_a3
integer :: s_i, i_i
! W1(j,i,a3,a) += (    1.00000000) T2(a1,a0,a,a2) D3(j,a3,i,a1,a2,a0) 
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_j = 0, nir-1
do s_a3 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_j,s_i) == IEOR(s_a3,s_a) .and. & 
IEOR(s_a1,s_a0) == IEOR(s_a,s_a2) .and. &
IEOR(IEOR(s_j,s_a3),s_i) == IEOR(IEOR(s_a1,s_a2),s_a0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a1,a0,a,a2) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a1, i_a0) =  &
  T2_(s_a, s_a0, s_a1)%array(i_a, i_a0, i_a1)
end do
end do
end do
! Z2 <-- D3(j,a3,i,a1,a2,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_j, i_a3, i_i) =  &
  D3_(s_a0, s_a2, s_a1, s_i, s_a3, s_j)%array(i_a0, i_a2, i_a1, i_i, i_a3, i_j)
end do
end do
end do
end do
end do

! Z3 <-- W1(j,i,a3,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W1(j,i,a3,a)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W1_(s_a, s_a3, s_i, s_j)%array(i_a, i_a3, i_i, i_j) = &
    W1_(s_a, s_a3, s_i, s_j)%array(i_a, i_a3, i_i, i_j) &
  + Z3_(i_a, i_j, i_a3, i_i)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

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

end subroutine g_sigma_oovv_ooov_no1_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x1_type0_noeri &
  (sb, ib, W1, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: W1(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xaaav(sleft, W1, nir, nsym, psym) ! -> Xaaav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x1_type0_noeri &
  (sb, ib, Xaaav, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xaaav)

end subroutine g_if_sigma_oovv_ooov_no1_x1_type0_noeri



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
subroutine g_sigma_oovv_ooov_no1_x1_type0_noeri &
  (s_b, i_b, W1_, S2_, Fc1_, nir, nsym, psym, flops)

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
type(symblock4), intent(inout) :: W1_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_j, i_j, s_i, i_i, s_a, i_a
! S2(i,j,a,b) += (    1.00000000) Fc1(b,a3) W1(j,i,a3,a) 
do s_a3 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_b,s_a3) == 0 .and. &
IEOR(s_j,s_i) == IEOR(s_a3,s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3) > 0) then

! Z1 <-- W1(j,i,a3,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a, i_a3) =  &
  W1_(s_a, s_a3, s_i, s_j)%array(i_a, i_a3, i_i, i_j)
end do
end do
end do
end do
! Z2 <-- Fc1(b,a3) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3) =  &
  Fc1_(s_a3, s_b)%array(i_a3, i_b)
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     1,&
                     psym(I_LENGTH,I_O, s_a3),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a))

! S2(i,j,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) = &
    S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) &
  + Z3_(i_j, i_i, i_a)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) * &
                1 * &
                psym(I_LENGTH,I_O, s_a3) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no1_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no2_x0_type0_noeri &
  (sa, ia, T2, W2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), W2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xaaa(sleft, W2, nir, nsym, psym) ! -> Xaaa (allocate) 
call g_sigma_oovv_ooov_no2_x0_type0_noeri &
  (sa, ia, av2_i, Xaaa, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_oovv_ooov_no2_x0_type0_noeri



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
subroutine g_sigma_oovv_ooov_no2_x0_type0_noeri &
  (s_a, i_a, T2_, W2_, D2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a0, i_a0, s_a2, i_a2, s_j, i_j, s_i, i_i
! W2(j,i,a2,a) += (    1.00000000) T2(a1,a0,a2,a) D2(j,a1,i,a0) 
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_a2 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_j,s_i) == IEOR(s_a2,s_a) .and. & 
IEOR(s_a1,s_a0) == IEOR(s_a2,s_a) .and. &
IEOR(s_j,s_a1) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a1,a0,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_a1, i_a0) =  &
  T2_(s_a2, s_a0, s_a1)%array(i_a2, i_a0, i_a1)
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

! Z3 <-- W2(j,i,a2,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a2))

! W2(j,i,a2,a)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W2_(s_a2, s_i, s_j)%array(i_a2, i_i, i_j) = &
    W2_(s_a2, s_i, s_j)%array(i_a2, i_i, i_j) &
  + Z3_(i_a2, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a2) * &
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

end subroutine g_sigma_oovv_ooov_no2_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no2_x1_type0_noeri &
  (sa, ia, sb, ib, W2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sb, ib
real(kind=8), intent(inout) :: W2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa

call set_symblock_Xaaa(sleft, W2, nir, nsym, psym) ! -> Xaaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no2_x1_type0_noeri &
  (sa, ia, sb, ib, Xaaa, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xaaa)

end subroutine g_if_sigma_oovv_ooov_no2_x1_type0_noeri



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
subroutine g_sigma_oovv_ooov_no2_x1_type0_noeri &
  (s_a, i_a, s_b, i_b, W2_, S2_, Fc1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_j, i_j, s_i, i_i
! S2(i,j,a,b) += (    1.00000000) Fc1(b,a2) W2(j,i,a2,a) 
do s_a2 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_b,s_a2) == 0 .and. &
IEOR(s_j,s_i) == IEOR(s_a2,s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- W2(j,i,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a2) =  &
  W2_(s_a2, s_i, s_j)%array(i_a2, i_i, i_j)
end do
end do
end do
! Z2 <-- Fc1(b,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2) =  &
  Fc1_(s_a2, s_b)%array(i_a2, i_b)
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(i,j,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) = &
    S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) &
  + Z3_(i_j, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
                1 * &
                psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no2_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no3_x0_type0_noeri &
  (sb, ib, T2, W3, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: T2(*), W3(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xaaa(sleft, W3, nir, nsym, psym) ! -> Xaaa (allocate) 
call g_sigma_oovv_ooov_no3_x0_type0_noeri &
  (sb, ib, av2_i, Xaaa, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_oovv_ooov_no3_x0_type0_noeri



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
subroutine g_sigma_oovv_ooov_no3_x0_type0_noeri &
  (s_b, i_b, T2_, W3_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W3_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_a1, i_a1, s_a2, i_a2, s_j, i_j, s_i, i_i
! W3(j,i,a2,b) += (    1.00000000) T2(a0,a1,a2,b) D2(j,a1,i,a0) 
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a2 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_j,s_i) == IEOR(s_a2,s_b) .and. & 
IEOR(s_a0,s_a1) == IEOR(s_a2,s_b) .and. &
IEOR(s_j,s_a1) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- T2(a0,a1,a2,b) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_a0, i_a1) =  &
  T2_(s_a2, s_a1, s_a0)%array(i_a2, i_a1, i_a0)
end do
end do
end do
! Z2 <-- D2(j,a1,i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a1, i_j, i_i) =  &
  D2_(s_a0, s_i, s_a1, s_j)%array(i_a0, i_i, i_a1, i_j)
end do
end do
end do
end do

! Z3 <-- W3(j,i,a2,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a2))

! W3(j,i,a2,b)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W3_(s_a2, s_i, s_j)%array(i_a2, i_i, i_j) = &
    W3_(s_a2, s_i, s_j)%array(i_a2, i_i, i_j) &
  + Z3_(i_a2, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a2) * &
                psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_oovv_ooov_no3_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no3_x1_type0_noeri &
  (sb, ib, W3, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: W3(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sb

call set_symblock_Xaaa(sleft, W3, nir, nsym, psym) ! -> Xaaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no3_x1_type0_noeri &
  (sb, ib, Xaaa, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xaaa)

end subroutine g_if_sigma_oovv_ooov_no3_x1_type0_noeri



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
subroutine g_sigma_oovv_ooov_no3_x1_type0_noeri &
  (s_b, i_b, W3_, S2_, Fc1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W3_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a, i_a, s_j, i_j, s_i, i_i
! S2(i,j,a,b) += (    1.00000000) Fc1(a2,a) W3(j,i,a2,b) 
do s_a2 = 0, nir-1
do s_a = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_a2,s_a) == 0 .and. &
IEOR(s_j,s_i) == IEOR(s_a2,s_b)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- Fc1(a2,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a2) =  &
  Fc1_(s_a, s_a2)%array(i_a, i_a2)
end do
end do
! Z2 <-- W3(j,i,a2,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_j, i_i) =  &
  W3_(s_a2, s_i, s_j)%array(i_a2, i_i, i_j)
end do
end do
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,j,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) = &
    S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) &
  + Z3_(i_a, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no3_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x0_type1_eri_o &
  (sa2, ia2, V2, W4, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa2, ia2
real(kind=8), intent(inout) :: V2(*), W4(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa2, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa2

call set_symblock_Xaaaav(sleft, W4, nir, nsym, psym) ! -> Xaaaav (allocate) 
call g_sigma_oovv_ooov_no0_x0_type1_eri_o &
  (sa2, ia2, h2_i, Xaaaav, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaaaav)

end subroutine g_if_sigma_oovv_ooov_no0_x0_type1_eri_o



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
! IS IT OK??
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Information for debugging 
! >> Score :: 111010
! RDM is rotated :: D3(j,a1,i,a3,a4,a0)  >> D3(a1,j,a3,i,a0,a4) 
! summedInd : @[a4, "active"] @[a3, "active"] 
! colInd : @[i, "active"] @[a0, "active"] @[j, "active"] @[a1, "active"] 
subroutine g_sigma_oovv_ooov_no0_x0_type1_eri_o &
  (s_a2, i_a2, V2_, W4_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a2, s_a2
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: W4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:,:,:)
real*8, allocatable :: Z3_(:,:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a4, i_a4, s_a3, i_a3, s_a, i_a, s_j, i_j, s_a1, i_a1
integer :: s_i, i_i, s_a0, i_a0
! W4(j,i,a1,a0,a2,a) += (    1.00000000) V2(a2,a4,a3,a) D3(j,a1,i,a3,a4,a0) 
do s_a4 = 0, nir-1
do s_a3 = 0, nir-1
do s_a = 0, nir-1
do s_j = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(IEOR(s_j,s_i),s_a1) == IEOR(IEOR(s_a0,s_a2),s_a) .and. & 
IEOR(s_a2,s_a4) == IEOR(s_a3,s_a) .and. &
IEOR(IEOR(s_j,s_a1),s_i) == IEOR(IEOR(s_a3,s_a4),s_a0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a3) > 0) then

! Z1 <-- V2(a2,a4,a3,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a4):psym(I_END,I_O, s_a4), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a4 = psym(I_BEGIN, I_O, s_a4), psym(I_END, I_O, s_a4)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a4, i_a3) =  &
  V2_(s_a, s_a3, s_a4)%array(i_a, i_a3, i_a4)
end do
end do
end do
! Z2 <-- D3(a1,j,a3,i,a0,a4) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a4):psym(I_END,I_O, s_a4), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a4 = psym(I_BEGIN, I_O, s_a4), psym(I_END, I_O, s_a4)
Z2_(i_a4, i_a3, i_i, i_a0, i_j, i_a1) =  &
  D3_(s_a4, s_a0, s_i, s_a3, s_j, s_a1)%array(i_a4, i_a0, i_i, i_a3, i_j, i_a1)
end do
end do
end do
end do
end do
end do

! Z3 <-- W4(j,i,a1,a0,a2,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a3),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W4(j,i,a1,a0,a2,a)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W4_(s_a, s_a0, s_a1, s_i, s_j)%array(i_a, i_a0, i_a1, i_i, i_j) = &
    W4_(s_a, s_a0, s_a1, s_i, s_j)%array(i_a, i_a0, i_a1, i_i, i_j) &
  + Z3_(i_a, i_i, i_a0, i_j, i_a1)
end do
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a3) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no0_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x1_type1_eri_o &
  (sa2, ia2, sb, ib, T2, W4, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa2, ia2, sb, ib
real(kind=8), intent(inout) :: T2(*), W4(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa2, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa2

call set_symblock_Xaaaav(sleft, W4, nir, nsym, psym) ! -> Xaaaav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no0_x1_type1_eri_o &
  (sa2, ia2, sb, ib, av2_i, Xaaaav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaaaav)

end subroutine g_if_sigma_oovv_ooov_no0_x1_type1_eri_o



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
subroutine g_sigma_oovv_ooov_no0_x1_type1_eri_o &
  (s_a2, i_a2, s_b, i_b, T2_, W4_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_a2, s_a2
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: W4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a0, i_a0, s_j, i_j, s_i, i_i, s_a, i_a
! S2(i,j,a,b) += (    1.00000000) T2(a1,a0,b,a2) W4(j,i,a1,a0,a2,a) 
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_a1,s_a0) == IEOR(s_b,s_a2) .and. &
IEOR(IEOR(s_j,s_i),s_a1) == IEOR(IEOR(s_a0,s_a2),s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W4(j,i,a1,a0,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a, i_a1, i_a0) =  &
  W4_(s_a, s_a0, s_a1, s_i, s_j)%array(i_a, i_a0, i_a1, i_i, i_j)
end do
end do
end do
end do
end do
! Z2 <-- T2(a1,a0,b,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  T2_(s_b, s_a0, s_a1)%array(i_b, i_a0, i_a1)
end do
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a))

! S2(i,j,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) = &
    S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) &
  + Z3_(i_j, i_i, i_a)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) * &
                1 * &
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

end subroutine g_sigma_oovv_ooov_no0_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x0_type1_eri_o &
  (sa3, ia3, V2, W7, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa3, ia3
real(kind=8), intent(inout) :: V2(*), W7(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa3, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa3

call set_symblock_Xaaaav(sleft, W7, nir, nsym, psym) ! -> Xaaaav (allocate) 
call g_sigma_oovv_ooov_no1_x0_type1_eri_o &
  (sa3, ia3, h2_i, Xaaaav, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaaaav)

end subroutine g_if_sigma_oovv_ooov_no1_x0_type1_eri_o



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
! IS IT OK??
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Information for debugging 
! >> Score :: 111111
! The best shot!
! RDM is rotated :: D3(j,a1,i,a0,a4,a2)  >> D3(a0,i,a1,j,a2,a4) 
! summedInd : @[a4, "active"] @[a2, "active"] 
! colInd : @[j, "active"] @[a1, "active"] @[i, "active"] @[a0, "active"] 
subroutine g_sigma_oovv_ooov_no1_x0_type1_eri_o &
  (s_a3, i_a3, V2_, W7_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a3, s_a3
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: W7_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_a4, i_a4, s_a2, i_a2, s_j, i_j, s_a1, i_a1
integer :: s_i, i_i, s_a0, i_a0
! W7(j,i,a1,a0,a3,a) += (    1.00000000) V2(a3,a,a4,a2) D3(j,a1,i,a0,a4,a2) 
do s_a = 0, nir-1
do s_a4 = 0, nir-1
do s_a2 = 0, nir-1
do s_j = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(IEOR(s_j,s_i),s_a1) == IEOR(IEOR(s_a0,s_a3),s_a) .and. & 
IEOR(s_a3,s_a) == IEOR(s_a4,s_a2) .and. &
IEOR(IEOR(s_j,s_a1),s_i) == IEOR(IEOR(s_a0,s_a4),s_a2)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(a3,a,a4,a2) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a4):psym(I_END,I_O, s_a4), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a4 = psym(I_BEGIN, I_O, s_a4), psym(I_END, I_O, s_a4)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a4, i_a2) =  &
  V2_(s_a2, s_a4, s_a)%array(i_a2, i_a4, i_a)
end do
end do
end do
! Z3 <-- W7(j,i,a1,a0,a3,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     D3_(s_a4, s_a2, s_j, s_a1, s_i, s_a0)%array,&
                     psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W7(j,i,a1,a0,a3,a)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W7_(s_a, s_a0, s_a1, s_i, s_j)%array(i_a, i_a0, i_a1, i_i, i_j) = &
    W7_(s_a, s_a0, s_a1, s_i, s_j)%array(i_a, i_a0, i_a1, i_i, i_j) &
  + Z3_(i_a, i_j, i_a1, i_i, i_a0)
end do
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no1_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x1_type1_eri_o &
  (sa3, ia3, sb, ib, T2, W7, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa3, ia3, sb, ib
real(kind=8), intent(inout) :: T2(*), W7(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa3, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa3

call set_symblock_Xaaaav(sleft, W7, nir, nsym, psym) ! -> Xaaaav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x1_type1_eri_o &
  (sa3, ia3, sb, ib, av2_i, Xaaaav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaaaav)

end subroutine g_if_sigma_oovv_ooov_no1_x1_type1_eri_o



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
subroutine g_sigma_oovv_ooov_no1_x1_type1_eri_o &
  (s_a3, i_a3, s_b, i_b, T2_, W7_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_a3, s_a3
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: W7_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a0, i_a0, s_j, i_j, s_i, i_i, s_a, i_a
! S2(i,j,a,b) += (    1.00000000) T2(a1,a0,b,a3) W7(j,i,a1,a0,a3,a) 
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_a1,s_a0) == IEOR(s_b,s_a3) .and. &
IEOR(IEOR(s_j,s_i),s_a1) == IEOR(IEOR(s_a0,s_a3),s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W7(j,i,a1,a0,a3,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a, i_a1, i_a0) =  &
  W7_(s_a, s_a0, s_a1, s_i, s_j)%array(i_a, i_a0, i_a1, i_i, i_j)
end do
end do
end do
end do
end do
! Z2 <-- T2(a1,a0,b,a3) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  T2_(s_b, s_a0, s_a1)%array(i_b, i_a0, i_a1)
end do
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a))

! S2(i,j,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) = &
    S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) &
  + Z3_(i_j, i_i, i_a)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) * &
                1 * &
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

end subroutine g_sigma_oovv_ooov_no1_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x0_type1_eri_v &
  (sa2, ia2, sb, ib, V2, W5, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa2, ia2, sb, ib
real(kind=8), intent(inout) :: V2(*), W5(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa2,sb)

call set_symblock_Xaaaa(sleft, W5, nir, nsym, psym) ! -> Xaaaa (allocate) 
call g_sigma_oovv_ooov_no0_x0_type1_eri_v &
  (sa2, ia2, sb, ib, h2_i, Xaaaa, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x0_type1_eri_v



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Information for debugging 
! >> Score :: 111010
! RDM is rotated :: D3(j,a3,i,a1,a4,a0)  >> D3(a0,a4,a3,j,a1,i) 
! rowInd : @[i, "active"] @[a1, "active"] @[j, "active"] @[a0, "active"] 
! summedInd : @[a4, "active"] @[a3, "active"] 
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_oovv_ooov_no0_x0_type1_eri_v &
  (s_a2, i_a2, s_b, i_b, V2_, W5_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_a2, s_a2
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W5_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_a4, i_a4, s_j, i_j, s_i, i_i, s_a1, i_a1
integer :: s_a0, i_a0
! W5(j,i,a1,a0,a2,b) += (    1.00000000) V2(b,a3,a4,a2) D3(j,a3,i,a1,a4,a0) 
do s_a3 = 0, nir-1
do s_a4 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(IEOR(s_j,s_i),s_a1) == IEOR(IEOR(s_a0,s_a2),s_b) .and. & 
IEOR(s_b,s_a3) == IEOR(s_a4,s_a2) .and. &
IEOR(IEOR(s_j,s_a3),s_i) == IEOR(IEOR(s_a1,s_a4),s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a3) > 0) then

! Z1 <-- D3(a0,a4,a3,j,a1,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a4):psym(I_END,I_O, s_a4), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a4 = psym(I_BEGIN, I_O, s_a4), psym(I_END, I_O, s_a4)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a1, i_j, i_a0, i_a4, i_a3) =  &
  D3_(s_i, s_a1, s_j, s_a3, s_a4, s_a0)%array(i_i, i_a1, i_j, i_a3, i_a4, i_a0)
end do
end do
end do
end do
end do
end do
! Z2 <-- V2(b,a3,a4,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a4):psym(I_END,I_O, s_a4), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a4 = psym(I_BEGIN, I_O, s_a4), psym(I_END, I_O, s_a4)
Z2_(i_a4, i_a3) =  &
  V2_(s_a2, s_a4, s_a3)%array(i_a2, i_a4, i_a3)
end do
end do

! Z3 <-- W5(j,i,a1,a0,a2,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a3),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a0))

! W5(j,i,a1,a0,a2,b)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W5_(s_a0, s_a1, s_i, s_j)%array(i_a0, i_a1, i_i, i_j) = &
    W5_(s_a0, s_a1, s_i, s_j)%array(i_a0, i_a1, i_i, i_j) &
  + Z3_(i_i, i_a1, i_j, i_a0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a3) * 2.0d+00

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

end subroutine g_sigma_oovv_ooov_no0_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x1_type1_eri_v &
  (sa2, ia2, sb, ib, T2, W5, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa2, ia2, sb, ib
real(kind=8), intent(inout) :: T2(*), W5(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa2, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa2,sb)

call set_symblock_Xaaaa(sleft, W5, nir, nsym, psym) ! -> Xaaaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no0_x1_type1_eri_v &
  (sa2, ia2, sb, ib, av2_i, Xaaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaaaa)

end subroutine g_if_sigma_oovv_ooov_no0_x1_type1_eri_v



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
subroutine g_sigma_oovv_ooov_no0_x1_type1_eri_v &
  (s_a2, i_a2, s_b, i_b, T2_, W5_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a2, s_a2
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W5_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a0, i_a0, s_a, i_a, s_j, i_j, s_i, i_i
! S2(i,j,a,b) += (    1.00000000) T2(a1,a0,a,a2) W5(j,i,a1,a0,a2,b) 
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_a1,s_a0) == IEOR(s_a,s_a2) .and. &
IEOR(IEOR(s_j,s_i),s_a1) == IEOR(IEOR(s_a0,s_a2),s_b)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a1,a0,a,a2) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a1, i_a0) =  &
  T2_(s_a, s_a0, s_a1)%array(i_a, i_a0, i_a1)
end do
end do
end do
! Z2 <-- W5(j,i,a1,a0,a2,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_j, i_i) =  &
  W5_(s_a0, s_a1, s_i, s_j)%array(i_a0, i_a1, i_i, i_j)
end do
end do
end do
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,j,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) = &
    S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) &
  + Z3_(i_a, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
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

end subroutine g_sigma_oovv_ooov_no0_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x0_type1_eri_v &
  (sb, ib, V2, W6, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: V2(*), W6(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sb

call set_symblock_Xaaaaa(sleft, W6, nir, nsym, psym) ! -> Xaaaaa (allocate) 
call g_sigma_oovv_ooov_no1_x0_type1_eri_v &
  (sb, ib, h2_i, Xaaaaa, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaaaaa)

end subroutine g_if_sigma_oovv_ooov_no1_x0_type1_eri_v



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
! IS IT OK??
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Information for debugging 
! >> Score :: 111111
! The best shot!
! RDM is rotated :: D3(j,a1,i,a0,a4,a2)  >> D3(a0,i,a1,j,a2,a4) 
! summedInd : @[a4, "active"] @[a2, "active"] 
! colInd : @[j, "active"] @[a1, "active"] @[i, "active"] @[a0, "active"] 
subroutine g_sigma_oovv_ooov_no1_x0_type1_eri_v &
  (s_b, i_b, V2_, W6_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: W6_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_a4, i_a4, s_a2, i_a2, s_j, i_j, s_a1, i_a1
integer :: s_i, i_i, s_a0, i_a0
! W6(j,i,a1,a0,a3,b) += (    1.00000000) V2(b,a3,a4,a2) D3(j,a1,i,a0,a4,a2) 
do s_a3 = 0, nir-1
do s_a4 = 0, nir-1
do s_a2 = 0, nir-1
do s_j = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(IEOR(s_j,s_i),s_a1) == IEOR(IEOR(s_a0,s_a3),s_b) .and. & 
IEOR(s_b,s_a3) == IEOR(s_a4,s_a2) .and. &
IEOR(IEOR(s_j,s_a1),s_i) == IEOR(IEOR(s_a0,s_a4),s_a2)) then

if(psym(I_LENGTH,I_O, s_a3) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(b,a3,a4,a2) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a4):psym(I_END,I_O, s_a4), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a4 = psym(I_BEGIN, I_O, s_a4), psym(I_END, I_O, s_a4)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z1_(i_a3, i_a4, i_a2) =  &
  V2_(s_a2, s_a4, s_a3)%array(i_a2, i_a4, i_a3)
end do
end do
end do
! Z3 <-- W6(j,i,a1,a0,a3,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a3),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a3),&
                     D3_(s_a4, s_a2, s_j, s_a1, s_i, s_a0)%array,&
                     psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a3))

! W6(j,i,a1,a0,a3,b)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
W6_(s_a3, s_a0, s_a1, s_i, s_j)%array(i_a3, i_a0, i_a1, i_i, i_j) = &
    W6_(s_a3, s_a0, s_a1, s_i, s_j)%array(i_a3, i_a0, i_a1, i_i, i_j) &
  + Z3_(i_a3, i_j, i_a1, i_i, i_a0)
end do
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a3) * &
                psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no1_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x1_type1_eri_v &
  (sa, ia, sb, ib, T2, W6, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sb, ib
real(kind=8), intent(inout) :: T2(*), W6(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xaaaaa(sleft, W6, nir, nsym, psym) ! -> Xaaaaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x1_type1_eri_v &
  (sa, ia, sb, ib, av2_i, Xaaaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaaaaa)

end subroutine g_if_sigma_oovv_ooov_no1_x1_type1_eri_v



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
subroutine g_sigma_oovv_ooov_no1_x1_type1_eri_v &
  (s_a, i_a, s_b, i_b, T2_, W6_, S2_, nir, nsym, psym, flops)

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
type(symblock5), intent(inout) :: W6_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a0, i_a0, s_a3, i_a3, s_j, i_j, s_i, i_i
! S2(i,j,a,b) += (    1.00000000) T2(a1,a0,a3,a) W6(j,i,a1,a0,a3,b) 
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_a3 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_a1,s_a0) == IEOR(s_a3,s_a) .and. &
IEOR(IEOR(s_j,s_i),s_a1) == IEOR(IEOR(s_a0,s_a3),s_b)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a3) > 0) then

! Z1 <-- W6(j,i,a1,a0,a3,b) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a1, i_a0, i_a3) =  &
  W6_(s_a3, s_a0, s_a1, s_i, s_j)%array(i_a3, i_a0, i_a1, i_i, i_j)
end do
end do
end do
end do
end do
! Z2 <-- T2(a1,a0,a3,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_a3) =  &
  T2_(s_a3, s_a0, s_a1)%array(i_a3, i_a0, i_a1)
end do
end do
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a3),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(i,j,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) = &
    S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) &
  + Z3_(i_j, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
                1 * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a3) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no1_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no2_x0_type1_eri_v &
  (sv0, iv0, T2, W8, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), W8(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xaaa(sleft, W8, nir, nsym, psym) ! -> Xaaa (allocate) 
call g_sigma_oovv_ooov_no2_x0_type1_eri_v &
  (sv0, iv0, av2_i, Xaaa, d3, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_oovv_ooov_no2_x0_type1_eri_v



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Information for debugging 
! >> Score :: 111111
! The best shot!
! RDM is rotated :: D3(j,a1,i,a3,a2,a0)  >> D3(a0,a2,a1,j,a3,i) 
! rowInd : @[i, "active"] @[a3, "active"] @[j, "active"] 
! summedInd : @[a1, "active"] @[a2, "active"] @[a0, "active"] 
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_oovv_ooov_no2_x0_type1_eri_v &
  (s_v0, i_v0, T2_, W8_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W8_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_a1, i_a1, s_a2, i_a2, s_j, i_j, s_i, i_i
integer :: s_a3, i_a3
! W8(j,i,a3,v0) += (    1.00000000) T2(a0,a1,a2,v0) D3(j,a1,i,a3,a2,a0) 
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a2 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a3 = 0, nir-1
if( &
IEOR(s_j,s_i) == IEOR(s_a3,s_v0) .and. & 
IEOR(s_a0,s_a1) == IEOR(s_a2,s_v0) .and. &
IEOR(IEOR(s_j,s_a1),s_i) == IEOR(IEOR(s_a3,s_a2),s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z2 <-- T2(a0,a1,a2,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a2, i_a0) =  &
  T2_(s_a2, s_a1, s_a0)%array(i_a2, i_a1, i_a0)
end do
end do
end do

! Z3 <-- W8(j,i,a3,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_j),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     D3_(s_i, s_a3, s_j, s_a1, s_a2, s_a0)%array,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_j))

! W8(j,i,a3,v0)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
W8_(s_a3, s_i, s_j)%array(i_a3, i_i, i_j) = &
    W8_(s_a3, s_i, s_j)%array(i_a3, i_i, i_j) &
  + Z3_(i_i, i_a3, i_j)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_j) * &
                1 * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

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

end subroutine g_sigma_oovv_ooov_no2_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no2_x1_type1_eri_v &
  (sb, ib, sv0, iv0, V2, W8, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W8(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xaaa(sleft, W8, nir, nsym, psym) ! -> Xaaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no2_x1_type1_eri_v &
  (sb, ib, sv0, iv0, h2_i, Xaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_oovv_ooov_no2_x1_type1_eri_v



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
subroutine g_sigma_oovv_ooov_no2_x1_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, V2_, W8_, S2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W8_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_a, i_a, s_j, i_j, s_i, i_i
! S2(i,j,a,b) += (    1.00000000) V2(v0,b,a3,a) W8(j,i,a3,v0) 
do s_a3 = 0, nir-1
do s_a = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_v0,s_b) == IEOR(s_a3,s_a) .and. &
IEOR(s_j,s_i) == IEOR(s_a3,s_v0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3) > 0) then

! Z1 <-- V2(v0,b,a3,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a3) =  &
  V2_(s_a, s_a3, s_b)%array(i_a, i_a3, i_b)
end do
end do
! Z2 <-- W8(j,i,a3,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_j, i_i) =  &
  W8_(s_a3, s_i, s_j)%array(i_a3, i_i, i_j)
end do
end do
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a3),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,j,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) = &
    S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) &
  + Z3_(i_a, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a3) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no2_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no3_x0_type1_eri_v &
  (sv0, iv0, T2, W9, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), W9(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xaaa(sleft, W9, nir, nsym, psym) ! -> Xaaa (allocate) 
call g_sigma_oovv_ooov_no3_x0_type1_eri_v &
  (sv0, iv0, av2_i, Xaaa, d3, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_oovv_ooov_no3_x0_type1_eri_v



!                  >> binary_contract3_new subroutine is called <<            
! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Information for debugging 
! >> Score :: 111111
! The best shot!
! RDM is rotated :: D3(j,a3,i,a1,a2,a0)  >> D3(a0,a2,a1,i,a3,j) 
! rowInd : @[j, "active"] @[a3, "active"] @[i, "active"] 
! summedInd : @[a1, "active"] @[a2, "active"] @[a0, "active"] 
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_oovv_ooov_no3_x0_type1_eri_v &
  (s_v0, i_v0, T2_, W9_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W9_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_a1, i_a1, s_a2, i_a2, s_j, i_j, s_a3, i_a3
integer :: s_i, i_i
! W9(j,i,a3,v0) += (    1.00000000) T2(a0,a1,a2,v0) D3(j,a3,i,a1,a2,a0) 
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a2 = 0, nir-1
do s_j = 0, nir-1
do s_a3 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_j,s_i) == IEOR(s_a3,s_v0) .and. & 
IEOR(s_a0,s_a1) == IEOR(s_a2,s_v0) .and. &
IEOR(IEOR(s_j,s_a3),s_i) == IEOR(IEOR(s_a1,s_a2),s_a0)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z2 <-- T2(a0,a1,a2,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a2, i_a0) =  &
  T2_(s_a2, s_a1, s_a0)%array(i_a2, i_a1, i_a0)
end do
end do
end do

! Z3 <-- W9(j,i,a3,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     D3_(s_j, s_a3, s_i, s_a1, s_a2, s_a0)%array,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_i))

! W9(j,i,a3,v0)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
W9_(s_a3, s_i, s_j)%array(i_a3, i_i, i_j) = &
    W9_(s_a3, s_i, s_j)%array(i_a3, i_i, i_j) &
  + Z3_(i_j, i_a3, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_i) * &
                1 * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

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

end subroutine g_sigma_oovv_ooov_no3_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no3_x1_type1_eri_v &
  (sb, ib, sv0, iv0, V2, W9, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W9(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xaaa(sleft, W9, nir, nsym, psym) ! -> Xaaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no3_x1_type1_eri_v &
  (sb, ib, sv0, iv0, h2_i, Xaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_oovv_ooov_no3_x1_type1_eri_v



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
subroutine g_sigma_oovv_ooov_no3_x1_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, V2_, W9_, S2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W9_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_a3, i_a3, s_j, i_j, s_i, i_i
! S2(i,j,a,b) += (    1.00000000) V2(v0,a,b,a3) W9(j,i,a3,v0) 
do s_a = 0, nir-1
do s_a3 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_v0,s_a) == IEOR(s_b,s_a3) .and. &
IEOR(s_j,s_i) == IEOR(s_a3,s_v0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3) > 0) then

! Z1 <-- V2(v0,a,b,a3) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a3) =  &
  V2_(s_a3, s_b, s_a)%array(i_a3, i_b, i_a)
end do
end do
! Z2 <-- W9(j,i,a3,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_j, i_i) =  &
  W9_(s_a3, s_i, s_j)%array(i_a3, i_i, i_j)
end do
end do
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a3),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,j,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) = &
    S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) &
  + Z3_(i_a, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a3) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no3_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no4_x0_type1_eri_v &
  (sa2, ia2, sb, ib, T2, V2, W10, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa2, ia2, sb, ib
real(kind=8), intent(inout) :: T2(*), V2(*), W10(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa2, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xaav(sleft, W10, nir, nsym, psym) ! -> Xaav (allocate) 
call g_sigma_oovv_ooov_no4_x0_type1_eri_v &
  (sa2, ia2, sb, ib, av2_i, h2_i, Xaav, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xaav)

end subroutine g_if_sigma_oovv_ooov_no4_x0_type1_eri_v



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
subroutine g_sigma_oovv_ooov_no4_x0_type1_eri_v &
  (s_a2, i_a2, s_b, i_b, T2_, V2_, W10_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_a2, s_a2
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W10_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_a, i_a, s_a1, i_a1, s_a0, i_a0
! W10(a1,a0,b,a) += (    1.00000000) V2(b,v0,a2,a) T2(a1,a0,v0,a2) 
do s_v0 = 0, nir-1
do s_a = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_a1,s_a0) == IEOR(s_b,s_a) .and. & 
IEOR(s_b,s_v0) == IEOR(s_a2,s_a) .and. &
IEOR(s_a1,s_a0) == IEOR(s_v0,s_a2)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- V2(b,v0,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_v0) =  &
  V2_(s_a, s_a2, s_v0)%array(i_a, i_a2, i_v0)
end do
end do
! Z2 <-- T2(a1,a0,v0,a2) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_a1, i_a0) =  &
  T2_(s_v0, s_a0, s_a1)%array(i_v0, i_a0, i_a1)
end do
end do
end do

! Z3 <-- W10(a1,a0,b,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W10(a1,a0,b,a)  <-- Z3
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W10_(s_a, s_a0, s_a1)%array(i_a, i_a0, i_a1) = &
    W10_(s_a, s_a0, s_a1)%array(i_a, i_a0, i_a1) &
  + Z3_(i_a, i_a1, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no4_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no4_x1_type1_eri_v &
  (sb, ib, W10, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: W10(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sb

call set_symblock_Xaav(sleft, W10, nir, nsym, psym) ! -> Xaav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no4_x1_type1_eri_v &
  (sb, ib, Xaav, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xaav)

end subroutine g_if_sigma_oovv_ooov_no4_x1_type1_eri_v



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
subroutine g_sigma_oovv_ooov_no4_x1_type1_eri_v &
  (s_b, i_b, W10_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W10_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_a1, i_a1, s_i, i_i, s_a0, i_a0, s_a, i_a
! S2(i,j,a,b) += (    1.00000000) D2(j,a1,i,a0) W10(a1,a0,b,a) 
do s_j = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_j,s_a1) == IEOR(s_i,s_a0) .and. &
IEOR(s_a1,s_a0) == IEOR(s_b,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W10(a1,a0,b,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a1, i_a0) =  &
  W10_(s_a, s_a0, s_a1)%array(i_a, i_a0, i_a1)
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

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,j,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) = &
    S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) &
  + Z3_(i_a, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
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

end subroutine g_sigma_oovv_ooov_no4_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no5_x0_type1_eri_v &
  (sb, ib, sv0, iv0, T2, V2, W11, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), W11(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xaav(sleft, W11, nir, nsym, psym) ! -> Xaav (allocate) 
call g_sigma_oovv_ooov_no5_x0_type1_eri_v &
  (sb, ib, sv0, iv0, av2_i, h2_i, Xaav, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xaav)

end subroutine g_if_sigma_oovv_ooov_no5_x0_type1_eri_v



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
subroutine g_sigma_oovv_ooov_no5_x0_type1_eri_v &
  (s_b, i_b, s_v0, i_v0, T2_, V2_, W11_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W11_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a, i_a, s_a1, i_a1, s_a0, i_a0
! W11(a1,a0,b,a) += (    1.00000000) V2(b,a2,v0,a) T2(a1,a0,a2,v0) 
do s_a2 = 0, nir-1
do s_a = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_a1,s_a0) == IEOR(s_b,s_a) .and. & 
IEOR(s_b,s_a2) == IEOR(s_v0,s_a) .and. &
IEOR(s_a1,s_a0) == IEOR(s_a2,s_v0)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(b,a2,v0,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a2) =  &
  V2_(s_a, s_v0, s_a2)%array(i_a, i_v0, i_a2)
end do
end do
! Z2 <-- T2(a1,a0,a2,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1, i_a0) =  &
  T2_(s_a2, s_a0, s_a1)%array(i_a2, i_a0, i_a1)
end do
end do
end do

! Z3 <-- W11(a1,a0,b,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W11(a1,a0,b,a)  <-- Z3
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W11_(s_a, s_a0, s_a1)%array(i_a, i_a0, i_a1) = &
    W11_(s_a, s_a0, s_a1)%array(i_a, i_a0, i_a1) &
  + Z3_(i_a, i_a1, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no5_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no5_x1_type1_eri_v &
  (sb, ib, W11, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: W11(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sb

call set_symblock_Xaav(sleft, W11, nir, nsym, psym) ! -> Xaav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no5_x1_type1_eri_v &
  (sb, ib, Xaav, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xaav)

end subroutine g_if_sigma_oovv_ooov_no5_x1_type1_eri_v



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
subroutine g_sigma_oovv_ooov_no5_x1_type1_eri_v &
  (s_b, i_b, W11_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W11_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_a1, i_a1, s_i, i_i, s_a0, i_a0, s_a, i_a
! S2(i,j,a,b) += (    1.00000000) D2(j,a1,i,a0) W11(a1,a0,b,a) 
do s_j = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_j,s_a1) == IEOR(s_i,s_a0) .and. &
IEOR(s_a1,s_a0) == IEOR(s_b,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W11(a1,a0,b,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a1, i_a0) =  &
  W11_(s_a, s_a0, s_a1)%array(i_a, i_a0, i_a1)
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

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,j,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) = &
    S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) &
  + Z3_(i_a, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
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

end subroutine g_sigma_oovv_ooov_no5_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no0_x0_type1_d4c_v &
  (sa, ia, sb, ib, C5, T2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sb, ib
real(kind=8), intent(inout) :: C5(*), T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_d4c(sa, C5, nir, nsym, psym) ! -> d4cf (allocate)
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no0_x0_type1_d4c_v &
  (sa, ia, sb, ib, d4cf, av2_i, av2_i2, nir, nsym, psym, flops)

deallocate(d4cf)
deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_oovv_ooov_no0_x0_type1_d4c_v



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
subroutine g_sigma_oovv_ooov_no0_x0_type1_d4c_v &
  (s_a, i_a, s_b, i_b, C5_, T2_, S2_, nir, nsym, psym, flops)

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
type(symblock5), intent(inout) :: C5_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_a1, i_a1, s_a2, i_a2, s_j, i_j, s_i, i_i
! S2(i,j,a,b) += (    1.00000000) T2(a0,a1,a2,b) C5(a2,a0,j,a1,i,a) 
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a2 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_a0,s_a1) == IEOR(s_a2,s_b) .and. &
IEOR(IEOR(s_a2,s_a0),s_j) == IEOR(IEOR(s_a1,s_i),s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- C5(a2,a0,j,a1,i,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a2, i_a0, i_a1) =  &
  C5_(s_i, s_a1, s_j, s_a0, s_a2)%array(i_i, i_a1, i_j, i_a0, i_a2)
end do
end do
end do
end do
end do
! Z2 <-- T2(a0,a1,a2,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a0, i_a1) =  &
  T2_(s_a2, s_a1, s_a0)%array(i_a2, i_a1, i_a0)
end do
end do
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(i,j,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) = &
    S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) &
  + Z3_(i_j, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
                1 * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_ooov_no0_x0_type1_d4c_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_ooov_no1_x0_type1_d4c_v &
  (sa2, ia2, sb, ib, C5, T2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa2, ia2, sb, ib
real(kind=8), intent(inout) :: C5(*), T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa2, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_d4c(sb, C5, nir, nsym, psym) ! -> d4cf (allocate)
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_ooov_no1_x0_type1_d4c_v &
  (sa2, ia2, sb, ib, d4cf, av2_i, av2_i2, nir, nsym, psym, flops)

deallocate(d4cf)
deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_oovv_ooov_no1_x0_type1_d4c_v



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
subroutine g_sigma_oovv_ooov_no1_x0_type1_d4c_v &
  (s_a2, i_a2, s_b, i_b, C5_, T2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a2, s_a2
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock5), intent(inout) :: C5_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a0, i_a0, s_a, i_a, s_i, i_i, s_j, i_j
! S2(i,j,a,b) += (    1.00000000) T2(a1,a0,a,a2) C5(a2,a0,i,a1,j,b) 
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
do s_j = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_a1,s_a0) == IEOR(s_a,s_a2) .and. &
IEOR(IEOR(s_a2,s_a0),s_i) == IEOR(IEOR(s_a1,s_j),s_b)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a1,a0,a,a2) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a1, i_a0) =  &
  T2_(s_a, s_a0, s_a1)%array(i_a, i_a0, i_a1)
end do
end do
end do
! Z2 <-- C5(a2,a0,i,a1,j,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j)))
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_i, i_j) =  &
  C5_(s_j, s_a1, s_i, s_a0, s_a2)%array(i_j, i_a1, i_i, i_a0, i_a2)
end do
end do
end do
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,j,a,b)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) = &
    S2_(s_a, s_j, s_i)%array(i_a, i_j, i_i) &
  + Z3_(i_a, i_i, i_j)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_oovv_ooov_no1_x0_type1_d4c_v

