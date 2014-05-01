#include <sci/icmr/fsrc/f_mr.fh>


!      ______                  __           
!     / ____/___   ____ ___   / /_ ____     
!    / /_   / _ \ / __ `__ \ / __// __ \ 
!   / __/  /  __// / / / / // /_ / /_/ /    
!  /_/     \___//_/ /_/ /_/ \__/ \____/  

!                                    Generated date : Sun Apr 20 10:26:10 2014



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_coov_no0_x0_type1_eri_c &
  (sc0, ic0, V2, W0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sc0, ic0
real(kind=8), intent(inout) :: V2(*), W0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sc0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc0

call set_symblock_Xaaaav(sleft, W0, nir, nsym, psym) ! -> Xaaaav (allocate) 
call g_sigma_oovv_coov_no0_x0_type1_eri_c &
  (sc0, ic0, h2_i, Xaaaav, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaaaav)

end subroutine g_if_sigma_oovv_coov_no0_x0_type1_eri_c



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
! >> Score :: 111001
! RDM is rotated :: D3(j,a0,i,a3,a1,a2)  >> D3(j,a0,i,a3,a1,a2) 
! summedInd : @[a2, "active"] @[a3, "active"] 
! colInd : @[a1, "active"] @[i, "active"] @[a0, "active"] @[j, "active"] 
subroutine g_sigma_oovv_coov_no0_x0_type1_eri_c &
  (s_c0, i_c0, V2_, W0_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c0, s_c0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:,:,:)
real*8, allocatable :: Z3_(:,:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a3, i_a3, s_a, i_a, s_j, i_j, s_a0, i_a0
integer :: s_i, i_i, s_a1, i_a1
! W0(c0,j,i,a1,a0,a) += (    1.00000000) V2(c0,a2,a3,a) D3(j,a0,i,a3,a1,a2) 
do s_a2 = 0, nir-1
do s_a3 = 0, nir-1
do s_a = 0, nir-1
do s_j = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(IEOR(s_c0,s_j),s_i) == IEOR(IEOR(s_a1,s_a0),s_a) .and. & 
IEOR(s_c0,s_a2) == IEOR(s_a3,s_a) .and. &
IEOR(IEOR(s_j,s_a0),s_i) == IEOR(IEOR(s_a3,s_a1),s_a2)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3) > 0) then

! Z1 <-- V2(c0,a2,a3,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a2, i_a3) =  &
  V2_(s_a, s_a3, s_a2)%array(i_a, i_a3, i_a2)
end do
end do
end do
! Z2 <-- D3(j,a0,i,a3,a1,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j)))
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a3, i_a1, i_i, i_a0, i_j) =  &
  D3_(s_a2, s_a1, s_a3, s_i, s_a0, s_j)%array(i_a2, i_a1, i_a3, i_i, i_a0, i_j)
end do
end do
end do
end do
end do
end do

! Z3 <-- W0(c0,j,i,a1,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W0(c0,j,i,a1,a0,a)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W0_(s_a, s_a0, s_a1, s_i, s_j)%array(i_a, i_a0, i_a1, i_i, i_j) = &
    W0_(s_a, s_a0, s_a1, s_i, s_j)%array(i_a, i_a0, i_a1, i_i, i_j) &
  + Z3_(i_a, i_a1, i_i, i_a0, i_j)
end do
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_j) * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3) * 2.0d+00

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

end subroutine g_sigma_oovv_coov_no0_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_coov_no0_x1_type1_eri_c &
  (sb, ib, sc0, ic0, T2, W0, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib, sc0, ic0
real(kind=8), intent(inout) :: T2(*), W0(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc0

call set_symblock_Xaaaav(sleft, W0, nir, nsym, psym) ! -> Xaaaav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_coov_no0_x1_type1_eri_c &
  (sb, ib, sc0, ic0, av2_i, Xaaaav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaaaav)

end subroutine g_if_sigma_oovv_coov_no0_x1_type1_eri_c



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
subroutine g_sigma_oovv_coov_no0_x1_type1_eri_c &
  (s_b, i_b, s_c0, i_c0, T2_, W0_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c0, s_c0
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_a1, i_a1, s_j, i_j, s_i, i_i, s_a, i_a
! S2(i,j,a,b) += (   -1.00000000) T2(c0,a0,a1,b) W0(c0,j,i,a1,a0,a) 
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_c0,s_a0) == IEOR(s_a1,s_b) .and. &
IEOR(IEOR(s_c0,s_j),s_i) == IEOR(IEOR(s_a1,s_a0),s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W0(c0,j,i,a1,a0,a) 
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
  W0_(s_a, s_a0, s_a1, s_i, s_j)%array(i_a, i_a0, i_a1, i_i, i_j)
end do
end do
end do
end do
end do
! Z2 <-- T2(c0,a0,a1,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  T2_(s_a1, s_a0, s_c0)%array(i_a1, i_a0, i_c0)
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
                     - 1.00000000d+00, &
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

end subroutine g_sigma_oovv_coov_no0_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_coov_no0_x0_type1_eri_o &
  (sa2, ia2, sb, ib, T2, V2, W2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa2, ia2, sb, ib
real(kind=8), intent(inout) :: T2(*), V2(*), W2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa2, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa2,sb)

call set_symblock_Xav(sleft, W2, nir, nsym, psym) ! -> Xav (allocate) 
call g_sigma_oovv_coov_no0_x0_type1_eri_o &
  (sa2, ia2, sb, ib, av2_i, h2_i, Xav, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xav)

end subroutine g_if_sigma_oovv_coov_no0_x0_type1_eri_o



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
subroutine g_sigma_oovv_coov_no0_x0_type1_eri_o &
  (s_a2, i_a2, s_b, i_b, T2_, V2_, W2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W2_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_c0, i_c0, s_a1, i_a1, s_a0, i_a0
! W2(a0,a2,b,a) += (    1.00000000) V2(a2,a,c0,a1) T2(c0,a0,a1,b) 
do s_a = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_a0,s_a2) == IEOR(s_b,s_a) .and. & 
IEOR(s_a2,s_a) == IEOR(s_c0,s_a1) .and. &
IEOR(s_c0,s_a0) == IEOR(s_a1,s_b)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(a2,a,c0,a1) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_c0, i_a1) =  &
  V2_(s_a1, s_c0, s_a)%array(i_a1, i_c0, i_a)
end do
end do
end do
! Z2 <-- T2(c0,a0,a1,b) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a1, i_a0) =  &
  T2_(s_a1, s_a0, s_c0)%array(i_a1, i_a0, i_c0)
end do
end do
end do

! Z3 <-- W2(a0,a2,b,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W2(a0,a2,b,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W2_(s_a, s_a0)%array(i_a, i_a0) = &
    W2_(s_a, s_a0)%array(i_a, i_a0) &
  + Z3_(i_a, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_coov_no0_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_coov_no0_x1_type1_eri_o &
  (sa2, ia2, sb, ib, W2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa2, ia2, sb, ib
real(kind=8), intent(inout) :: W2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa2,sb)

call set_symblock_Xav(sleft, W2, nir, nsym, psym) ! -> Xav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_coov_no0_x1_type1_eri_o &
  (sa2, ia2, sb, ib, Xav, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xav)

end subroutine g_if_sigma_oovv_coov_no0_x1_type1_eri_o



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
subroutine g_sigma_oovv_coov_no0_x1_type1_eri_o &
  (s_a2, i_a2, s_b, i_b, W2_, S2_, D2_, nir, nsym, psym, flops)

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
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W2_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_a0, i_a0, s_i, i_i, s_a, i_a
! S2(i,j,a,b) += (    2.00000000) D2(j,a0,i,a2) W2(a0,a2,b,a) 
do s_j = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_j,s_a0) == IEOR(s_i,s_a2) .and. &
IEOR(s_a0,s_a2) == IEOR(s_b,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W2(a0,a2,b,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0) =  &
  W2_(s_a, s_a0)%array(i_a, i_a0)
end do
end do
! Z2 <-- D2(j,a0,i,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_j, i_i) =  &
  D2_(s_a2, s_i, s_a0, s_j)%array(i_a2, i_i, i_a0, i_j)
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
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
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
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_coov_no0_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_coov_no1_x0_type1_eri_o &
  (sa1, ia1, sb, ib, T2, V2, W5, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa1, ia1, sb, ib
real(kind=8), intent(inout) :: T2(*), V2(*), W5(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa1,sb)

call set_symblock_Xav(sleft, W5, nir, nsym, psym) ! -> Xav (allocate) 
call g_sigma_oovv_coov_no1_x0_type1_eri_o &
  (sa1, ia1, sb, ib, av2_i, h2_i, Xav, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xav)

end subroutine g_if_sigma_oovv_coov_no1_x0_type1_eri_o



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
subroutine g_sigma_oovv_coov_no1_x0_type1_eri_o &
  (s_a1, i_a1, s_b, i_b, T2_, V2_, W5_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W5_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a2, i_a2, s_a, i_a, s_a0, i_a0
! W5(a0,a1,b,a) += (    1.00000000) V2(a1,c0,a2,a) T2(c0,a0,a2,b) 
do s_c0 = 0, nir-1
do s_a2 = 0, nir-1
do s_a = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_a0,s_a1) == IEOR(s_b,s_a) .and. & 
IEOR(s_a1,s_c0) == IEOR(s_a2,s_a) .and. &
IEOR(s_c0,s_a0) == IEOR(s_a2,s_b)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(a1,c0,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_c0, i_a2) =  &
  V2_(s_a, s_a2, s_c0)%array(i_a, i_a2, i_c0)
end do
end do
end do
! Z2 <-- T2(c0,a0,a2,b) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a2, i_a0) =  &
  T2_(s_a2, s_a0, s_c0)%array(i_a2, i_a0, i_c0)
end do
end do
end do

! Z3 <-- W5(a0,a1,b,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W5(a0,a1,b,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
W5_(s_a, s_a0)%array(i_a, i_a0) = &
    W5_(s_a, s_a0)%array(i_a, i_a0) &
  + Z3_(i_a, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_coov_no1_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_coov_no1_x1_type1_eri_o &
  (sa1, ia1, sb, ib, W5, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa1, ia1, sb, ib
real(kind=8), intent(inout) :: W5(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa1,sb)

call set_symblock_Xav(sleft, W5, nir, nsym, psym) ! -> Xav (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_coov_no1_x1_type1_eri_o &
  (sa1, ia1, sb, ib, Xav, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xav)

end subroutine g_if_sigma_oovv_coov_no1_x1_type1_eri_o



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
subroutine g_sigma_oovv_coov_no1_x1_type1_eri_o &
  (s_a1, i_a1, s_b, i_b, W5_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W5_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_a0, i_a0, s_i, i_i, s_a, i_a
! S2(i,j,a,b) += (   -1.00000000) D2(j,a0,i,a1) W5(a0,a1,b,a) 
do s_j = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_j,s_a0) == IEOR(s_i,s_a1) .and. &
IEOR(s_a0,s_a1) == IEOR(s_b,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W5(a0,a1,b,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_a, i_a0) =  &
  W5_(s_a, s_a0)%array(i_a, i_a0)
end do
end do
! Z2 <-- D2(j,a0,i,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_j, i_i) =  &
  D2_(s_a1, s_i, s_a0, s_j)%array(i_a1, i_i, i_a0, i_j)
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
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
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
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_coov_no1_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_coov_no0_x0_type1_eri_v &
  (sb, ib, V2, W1, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: V2(*), W1(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sb

call set_symblock_Xcaaaa(sleft, W1, nir, nsym, psym) ! -> Xcaaaa (allocate) 
call g_sigma_oovv_coov_no0_x0_type1_eri_v &
  (sb, ib, h2_i, Xcaaaa, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaaaa)

end subroutine g_if_sigma_oovv_coov_no0_x0_type1_eri_v



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
! >> Score :: 110110
! RDM is rotated :: D3(j,a2,i,a0,a1,a3)  >> D3(a3,a1,a2,j,a0,i) 
! rowInd : @[i, "active"] @[a0, "active"] @[j, "active"] @[a1, "active"] 
! summedInd : @[a2, "active"] @[a3, "active"] 
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_oovv_coov_no0_x0_type1_eri_v &
  (s_b, i_b, V2_, W1_, D3_, nir, nsym, psym, flops)

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
type(symblock5), intent(inout) :: W1_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_c0, i_c0, s_a3, i_a3, s_j, i_j, s_i, i_i
integer :: s_a0, i_a0, s_a1, i_a1
! W1(c0,j,i,a1,a0,b) += (    1.00000000) V2(b,a2,c0,a3) D3(j,a2,i,a0,a1,a3) 
do s_a2 = 0, nir-1
do s_c0 = 0, nir-1
do s_a3 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(IEOR(s_c0,s_j),s_i) == IEOR(IEOR(s_a1,s_a0),s_b) .and. & 
IEOR(s_b,s_a2) == IEOR(s_c0,s_a3) .and. &
IEOR(IEOR(s_j,s_a2),s_i) == IEOR(IEOR(s_a0,s_a1),s_a3)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3) > 0) then

! Z1 <-- D3(a3,a1,a2,j,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_j, i_a1, i_a2, i_a3) =  &
  D3_(s_i, s_a0, s_j, s_a2, s_a1, s_a3)%array(i_i, i_a0, i_j, i_a2, i_a1, i_a3)
end do
end do
end do
end do
end do
end do
! Z2 <-- V2(b,a2,c0,a3) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a3, i_c0) =  &
  V2_(s_a3, s_c0, s_a2)%array(i_a3, i_c0, i_a2)
end do
end do
end do

! Z3 <-- W1(c0,j,i,a1,a0,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a1))

! W1(c0,j,i,a1,a0,b)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W1_(s_a0, s_a1, s_i, s_j, s_c0)%array(i_a0, i_a1, i_i, i_j, i_c0) = &
    W1_(s_a0, s_a1, s_i, s_j, s_c0)%array(i_a0, i_a1, i_i, i_j, i_c0) &
  + Z3_(i_i, i_a0, i_j, i_a1, i_c0)
end do
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3) * 2.0d+00

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

end subroutine g_sigma_oovv_coov_no0_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_coov_no0_x1_type1_eri_v &
  (sa, ia, sb, ib, T2, W1, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sb, ib
real(kind=8), intent(inout) :: T2(*), W1(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_Xcaaaa(sleft, W1, nir, nsym, psym) ! -> Xcaaaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_coov_no0_x1_type1_eri_v &
  (sa, ia, sb, ib, av2_i, Xcaaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaaaa)

end subroutine g_if_sigma_oovv_coov_no0_x1_type1_eri_v



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
subroutine g_sigma_oovv_coov_no0_x1_type1_eri_v &
  (s_a, i_a, s_b, i_b, T2_, W1_, S2_, nir, nsym, psym, flops)

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
type(symblock5), intent(inout) :: W1_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a0, i_a0, s_a1, i_a1, s_j, i_j, s_i, i_i
! S2(i,j,a,b) += (   -1.00000000) T2(c0,a0,a1,a) W1(c0,j,i,a1,a0,b) 
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_c0,s_a0) == IEOR(s_a1,s_a) .and. &
IEOR(IEOR(s_c0,s_j),s_i) == IEOR(IEOR(s_a1,s_a0),s_b)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W1(c0,j,i,a1,a0,b) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_c0, i_a1, i_a0) =  &
  W1_(s_a0, s_a1, s_i, s_j, s_c0)%array(i_a0, i_a1, i_i, i_j, i_c0)
end do
end do
end do
end do
end do
! Z2 <-- T2(c0,a0,a1,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a1, i_a0) =  &
  T2_(s_a1, s_a0, s_c0)%array(i_a1, i_a0, i_c0)
end do
end do
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
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
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_coov_no0_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_coov_no1_x0_type1_eri_v &
  (sa, ia, sb, ib, T2, V2, W3, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sb, ib
real(kind=8), intent(inout) :: T2(*), V2(*), W3(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa,sb)

call set_symblock_Xaa(sleft, W3, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_oovv_coov_no1_x0_type1_eri_v &
  (sa, ia, sb, ib, av2_i, h2_i, Xaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_oovv_coov_no1_x0_type1_eri_v



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
subroutine g_sigma_oovv_coov_no1_x0_type1_eri_v &
  (s_a, i_a, s_b, i_b, T2_, V2_, W3_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W3_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_c0, i_c0, s_a2, i_a2, s_a0, i_a0
! W3(a0,a1,a,b) += (    1.00000000) V2(b,a1,c0,a2) T2(c0,a0,a2,a) 
do s_a1 = 0, nir-1
do s_c0 = 0, nir-1
do s_a2 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_a0,s_a1) == IEOR(s_a,s_b) .and. & 
IEOR(s_b,s_a1) == IEOR(s_c0,s_a2) .and. &
IEOR(s_c0,s_a0) == IEOR(s_a2,s_a)) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(b,a1,c0,a2) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_c0, i_a2) =  &
  V2_(s_a2, s_c0, s_a1)%array(i_a2, i_c0, i_a1)
end do
end do
end do
! Z2 <-- T2(c0,a0,a2,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a2, i_a0) =  &
  T2_(s_a2, s_a0, s_c0)%array(i_a2, i_a0, i_c0)
end do
end do
end do

! Z3 <-- W3(a0,a1,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W3(a0,a1,a,b)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W3_(s_a1, s_a0)%array(i_a1, i_a0) = &
    W3_(s_a1, s_a0)%array(i_a1, i_a0) &
  + Z3_(i_a1, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_coov_no1_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_coov_no1_x1_type1_eri_v &
  (sa, ia, sb, ib, W3, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sb, ib
real(kind=8), intent(inout) :: W3(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa,sb)

call set_symblock_Xaa(sleft, W3, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_coov_no1_x1_type1_eri_v &
  (sa, ia, sb, ib, Xaa, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xaa)

end subroutine g_if_sigma_oovv_coov_no1_x1_type1_eri_v



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
subroutine g_sigma_oovv_coov_no1_x1_type1_eri_v &
  (s_a, i_a, s_b, i_b, W3_, S2_, D2_, nir, nsym, psym, flops)

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
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W3_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! S2(i,j,a,b) += (    2.00000000) D2(j,a1,i,a0) W3(a0,a1,a,b) 
do s_j = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_j,s_a1) == IEOR(s_i,s_a0) .and. &
IEOR(s_a0,s_a1) == IEOR(s_a,s_b)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(j,a1,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a1, i_a0) =  &
  D2_(s_a0, s_i, s_a1, s_j)%array(i_a0, i_i, i_a1, i_j)
end do
end do
end do
end do
! Z2 <-- W3(a0,a1,a,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  W3_(s_a1, s_a0)%array(i_a1, i_a0)
end do
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
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
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_coov_no1_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_coov_no2_x0_type1_eri_v &
  (sa, ia, sb, ib, T2, V2, W4, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sb, ib
real(kind=8), intent(inout) :: T2(*), V2(*), W4(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa,sb)

call set_symblock_Xaa(sleft, W4, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_oovv_coov_no2_x0_type1_eri_v &
  (sa, ia, sb, ib, av2_i, h2_i, Xaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_oovv_coov_no2_x0_type1_eri_v



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
subroutine g_sigma_oovv_coov_no2_x0_type1_eri_v &
  (s_a, i_a, s_b, i_b, T2_, V2_, W4_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W4_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_c0, i_c0, s_a2, i_a2, s_a0, i_a0
! W4(a0,a2,a,b) += (    1.00000000) V2(b,a1,c0,a2) T2(c0,a0,a1,a) 
do s_a1 = 0, nir-1
do s_c0 = 0, nir-1
do s_a2 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_a0,s_a2) == IEOR(s_a,s_b) .and. & 
IEOR(s_b,s_a1) == IEOR(s_c0,s_a2) .and. &
IEOR(s_c0,s_a0) == IEOR(s_a1,s_a)) then

if(psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- V2(b,a1,c0,a2) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_a1, i_c0) =  &
  V2_(s_a2, s_c0, s_a1)%array(i_a2, i_c0, i_a1)
end do
end do
end do
! Z2 <-- T2(c0,a0,a1,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_c0, i_a0) =  &
  T2_(s_a1, s_a0, s_c0)%array(i_a1, i_a0, i_c0)
end do
end do
end do

! Z3 <-- W4(a0,a2,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a2))

! W4(a0,a2,a,b)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W4_(s_a2, s_a0)%array(i_a2, i_a0) = &
    W4_(s_a2, s_a0)%array(i_a2, i_a0) &
  + Z3_(i_a2, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a2) * &
                psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_coov_no2_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_coov_no2_x1_type1_eri_v &
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

sleft = IEOR(sa,sb)

call set_symblock_Xaa(sleft, W4, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_coov_no2_x1_type1_eri_v &
  (sa, ia, sb, ib, Xaa, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xaa)

end subroutine g_if_sigma_oovv_coov_no2_x1_type1_eri_v



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
subroutine g_sigma_oovv_coov_no2_x1_type1_eri_v &
  (s_a, i_a, s_b, i_b, W4_, S2_, D2_, nir, nsym, psym, flops)

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
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W4_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_a2, i_a2, s_i, i_i, s_a0, i_a0
! S2(i,j,a,b) += (   -1.00000000) D2(j,a2,i,a0) W4(a0,a2,a,b) 
do s_j = 0, nir-1
do s_a2 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_a,s_b) .and. & 
IEOR(s_j,s_a2) == IEOR(s_i,s_a0) .and. &
IEOR(s_a0,s_a2) == IEOR(s_a,s_b)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(j,a2,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a2, i_a0) =  &
  D2_(s_a0, s_i, s_a2, s_j)%array(i_a0, i_i, i_a2, i_j)
end do
end do
end do
end do
! Z2 <-- W4(a0,a2,a,b) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a0) =  &
  W4_(s_a2, s_a0)%array(i_a2, i_a0)
end do
end do

! Z3 <-- S2(i,j,a,b) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0),&
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
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_coov_no2_x1_type1_eri_v

