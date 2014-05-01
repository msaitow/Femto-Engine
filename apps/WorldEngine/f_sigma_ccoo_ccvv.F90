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
subroutine g_if_sigma_ccoo_ccvv_no0_x0_type1_eri_v &
  (sj, ij, sv1, iv1, V2, W0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv1, iv1
real(kind=8), intent(inout) :: V2(*), W0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sj,sv1)

call set_symblock_Xav(sleft, W0, nir, nsym, psym) ! -> Xav (allocate) 
call g_sigma_ccoo_ccvv_no0_x0_type1_eri_v &
  (sj, ij, sv1, iv1, h2_i, Xav, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xav)

end subroutine g_if_sigma_ccoo_ccvv_no0_x0_type1_eri_v



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
subroutine g_sigma_ccoo_ccvv_no0_x0_type1_eri_v &
  (s_j, i_j, s_v1, i_v1, V2_, W0_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W0_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_v0, i_v0, s_a0, i_a0, s_i, i_i
! W0(j,i,v1,v0) += (    1.00000000) V2(v1,a1,v0,a0) D2(j,a1,i,a0) 
do s_a1 = 0, nir-1
do s_v0 = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_j,s_i) == IEOR(s_v1,s_v0) .and. & 
IEOR(s_v1,s_a1) == IEOR(s_v0,s_a0) .and. &
IEOR(s_j,s_a1) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v1,a1,v0,a0) 
allocate(Z1_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z1_(i_v0, i_a1, i_a0) =  &
  V2_(s_a0, s_v0, s_a1)%array(i_a0, i_v0, i_a1)
end do
end do
end do
! Z2 <-- D2(j,a1,i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_i) =  &
  D2_(s_a0, s_i, s_a1, s_j)%array(i_a0, i_i, i_a1, i_j)
end do
end do
end do

! Z3 <-- W0(j,i,v1,v0) 
allocate(Z3_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_v0),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_v0))

! W0(j,i,v1,v0)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
W0_(s_v0, s_i)%array(i_v0, i_i) = &
    W0_(s_v0, s_i)%array(i_v0, i_i) &
  + Z3_(i_v0, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_v0) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccoo_ccvv_no0_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccoo_ccvv_no0_x1_type1_eri_v &
  (sj, ij, sv1, iv1, T2, W0, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv1, iv1
real(kind=8), intent(inout) :: T2(*), W0(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sv1)

call set_symblock_Xav(sleft, W0, nir, nsym, psym) ! -> Xav (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccoo_ccvv_no0_x1_type1_eri_v &
  (sj, ij, sv1, iv1, av2_i, Xav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xav)

end subroutine g_if_sigma_ccoo_ccvv_no0_x1_type1_eri_v



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
subroutine g_sigma_ccoo_ccvv_no0_x1_type1_eri_v &
  (s_j, i_j, s_v1, i_v1, T2_, W0_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W0_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_x, i_x, s_v0, i_v0, s_i, i_i
! S2(w,x,i,j) += (    2.00000000) T2(w,x,v0,v1) W0(j,i,v1,v0) 
do s_w = 0, nir-1
do s_x = 0, nir-1
do s_v0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_x) == IEOR(s_v0,s_v1) .and. &
IEOR(s_j,s_i) == IEOR(s_v1,s_v0)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W0(j,i,v1,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_v0) =  &
  W0_(s_v0, s_i)%array(i_v0, i_i)
end do
end do
! Z2 <-- T2(w,x,v0,v1) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_w, i_x) =  &
  T2_(s_v0, s_x, s_w)%array(i_v0, i_x, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_V, s_v0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,j)  <-- Z3
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
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccoo_ccvv_no0_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccoo_ccvv_no1_x0_type1_eri_v &
  (sj, ij, sv1, iv1, V2, W1, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv1, iv1
real(kind=8), intent(inout) :: V2(*), W1(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sj,sv1)

call set_symblock_Xav(sleft, W1, nir, nsym, psym) ! -> Xav (allocate) 
call g_sigma_ccoo_ccvv_no1_x0_type1_eri_v &
  (sj, ij, sv1, iv1, h2_i, Xav, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xav)

end subroutine g_if_sigma_ccoo_ccvv_no1_x0_type1_eri_v



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
subroutine g_sigma_ccoo_ccvv_no1_x0_type1_eri_v &
  (s_j, i_j, s_v1, i_v1, V2_, W1_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W1_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_v0, i_v0, s_a0, i_a0
! W1(j,i,v1,v0) += (    1.00000000) V2(v1,i,v0,a0) D1(j,a0) 
do s_i = 0, nir-1
do s_v0 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_j,s_i) == IEOR(s_v1,s_v0) .and. & 
IEOR(s_v1,s_i) == IEOR(s_v0,s_a0) .and. &
IEOR(s_j,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v1,i,v0,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_v0, i_a0) =  &
  V2_(s_a0, s_v0, s_i)%array(i_a0, i_v0, i_i)
end do
end do
end do
! Z2 <-- D1(j,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0) =  &
  D1_(s_a0, s_j)%array(i_a0, i_j)
end do

! Z3 <-- W1(j,i,v1,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_v0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_v0))

! W1(j,i,v1,v0)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
W1_(s_v0, s_i)%array(i_v0, i_i) = &
    W1_(s_v0, s_i)%array(i_v0, i_i) &
  + Z3_(i_i, i_v0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_v0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccoo_ccvv_no1_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccoo_ccvv_no1_x1_type1_eri_v &
  (sj, ij, sv1, iv1, T2, W1, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv1, iv1
real(kind=8), intent(inout) :: T2(*), W1(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sv1)

call set_symblock_Xav(sleft, W1, nir, nsym, psym) ! -> Xav (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccoo_ccvv_no1_x1_type1_eri_v &
  (sj, ij, sv1, iv1, av2_i, Xav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xav)

end subroutine g_if_sigma_ccoo_ccvv_no1_x1_type1_eri_v



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
subroutine g_sigma_ccoo_ccvv_no1_x1_type1_eri_v &
  (s_j, i_j, s_v1, i_v1, T2_, W1_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W1_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_x, i_x, s_v0, i_v0, s_i, i_i
! S2(w,x,i,j) += (    2.00000000) T2(w,x,v0,v1) W1(j,i,v1,v0) 
do s_w = 0, nir-1
do s_x = 0, nir-1
do s_v0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_x) == IEOR(s_v0,s_v1) .and. &
IEOR(s_j,s_i) == IEOR(s_v1,s_v0)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W1(j,i,v1,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_v0) =  &
  W1_(s_v0, s_i)%array(i_v0, i_i)
end do
end do
! Z2 <-- T2(w,x,v0,v1) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_w, i_x) =  &
  T2_(s_v0, s_x, s_w)%array(i_v0, i_x, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_V, s_v0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,j)  <-- Z3
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
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccoo_ccvv_no1_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccoo_ccvv_no2_x0_type1_eri_v &
  (sj, ij, sv1, iv1, V2, W2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv1, iv1
real(kind=8), intent(inout) :: V2(*), W2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sj,sv1)

call set_symblock_Xav(sleft, W2, nir, nsym, psym) ! -> Xav (allocate) 
call g_sigma_ccoo_ccvv_no2_x0_type1_eri_v &
  (sj, ij, sv1, iv1, h2_i, Xav, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xav)

end subroutine g_if_sigma_ccoo_ccvv_no2_x0_type1_eri_v



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
subroutine g_sigma_ccoo_ccvv_no2_x0_type1_eri_v &
  (s_j, i_j, s_v1, i_v1, V2_, W2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W2_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_v0, i_v0, s_a0, i_a0
! W2(j,i,v1,v0) += (    1.00000000) V2(v1,i,v0,a0) D1(j,a0) 
do s_i = 0, nir-1
do s_v0 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_j,s_i) == IEOR(s_v1,s_v0) .and. & 
IEOR(s_v1,s_i) == IEOR(s_v0,s_a0) .and. &
IEOR(s_j,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v1,i,v0,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_v0, i_a0) =  &
  V2_(s_a0, s_v0, s_i)%array(i_a0, i_v0, i_i)
end do
end do
end do
! Z2 <-- D1(j,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0) =  &
  D1_(s_a0, s_j)%array(i_a0, i_j)
end do

! Z3 <-- W2(j,i,v1,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_v0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_v0))

! W2(j,i,v1,v0)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
W2_(s_v0, s_i)%array(i_v0, i_i) = &
    W2_(s_v0, s_i)%array(i_v0, i_i) &
  + Z3_(i_i, i_v0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_V, s_v0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccoo_ccvv_no2_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccoo_ccvv_no2_x1_type1_eri_v &
  (sj, ij, sv1, iv1, T2, W2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv1, iv1
real(kind=8), intent(inout) :: T2(*), W2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sv1)

call set_symblock_Xav(sleft, W2, nir, nsym, psym) ! -> Xav (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccoo_ccvv_no2_x1_type1_eri_v &
  (sj, ij, sv1, iv1, av2_i, Xav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xav)

end subroutine g_if_sigma_ccoo_ccvv_no2_x1_type1_eri_v



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
subroutine g_sigma_ccoo_ccvv_no2_x1_type1_eri_v &
  (s_j, i_j, s_v1, i_v1, T2_, W2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W2_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_v0, i_v0, s_i, i_i
! S2(w,x,i,j) += (   -4.00000000) T2(x,w,v0,v1) W2(j,i,v1,v0) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_v0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_j) .and. & 
IEOR(s_x,s_w) == IEOR(s_v0,s_v1) .and. &
IEOR(s_j,s_i) == IEOR(s_v1,s_v0)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W2(j,i,v1,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_v0) =  &
  W2_(s_v0, s_i)%array(i_v0, i_i)
end do
end do
! Z2 <-- T2(x,w,v0,v1) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_x, i_w) =  &
  T2_(s_v0, s_w, s_x)%array(i_v0, i_w, i_x)
end do
end do
end do

! Z3 <-- S2(w,x,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_V, s_v0),&
                     - 4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,j)  <-- Z3
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
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccoo_ccvv_no2_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccoo_ccvv_no3_x0_type1_eri_v &
  (sj, ij, sv1, iv1, V2, W3, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv1, iv1
real(kind=8), intent(inout) :: V2(*), W3(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sj,sv1)

call set_symblock_Xav(sleft, W3, nir, nsym, psym) ! -> Xav (allocate) 
call g_sigma_ccoo_ccvv_no3_x0_type1_eri_v &
  (sj, ij, sv1, iv1, h2_i, Xav, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xav)

end subroutine g_if_sigma_ccoo_ccvv_no3_x0_type1_eri_v



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
subroutine g_sigma_ccoo_ccvv_no3_x0_type1_eri_v &
  (s_j, i_j, s_v1, i_v1, V2_, W3_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W3_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_a0, i_a0, s_i, i_i
! W3(i,j,v1,v0) += (    1.00000000) V2(v1,j,v0,a0) D1(i,a0) 
do s_v0 = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_v1,s_v0) .and. & 
IEOR(s_v1,s_j) == IEOR(s_v0,s_a0) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v1,j,v0,a0) 
allocate(Z1_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z1_(i_v0, i_a0) =  &
  V2_(s_a0, s_v0, s_j)%array(i_a0, i_v0, i_j)
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

! Z3 <-- W3(i,j,v1,v0) 
allocate(Z3_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_v0),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_v0))

! W3(i,j,v1,v0)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
W3_(s_v0, s_i)%array(i_v0, i_i) = &
    W3_(s_v0, s_i)%array(i_v0, i_i) &
  + Z3_(i_v0, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_v0) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccoo_ccvv_no3_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccoo_ccvv_no3_x1_type1_eri_v &
  (sj, ij, sv1, iv1, T2, W3, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv1, iv1
real(kind=8), intent(inout) :: T2(*), W3(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sv1)

call set_symblock_Xav(sleft, W3, nir, nsym, psym) ! -> Xav (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccoo_ccvv_no3_x1_type1_eri_v &
  (sj, ij, sv1, iv1, av2_i, Xav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xav)

end subroutine g_if_sigma_ccoo_ccvv_no3_x1_type1_eri_v



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
subroutine g_sigma_ccoo_ccvv_no3_x1_type1_eri_v &
  (s_j, i_j, s_v1, i_v1, T2_, W3_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W3_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_x, i_x, s_v0, i_v0, s_i, i_i
! S2(w,x,i,j) += (   -4.00000000) T2(w,x,v0,v1) W3(i,j,v1,v0) 
do s_w = 0, nir-1
do s_x = 0, nir-1
do s_v0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_j) .and. & 
IEOR(s_w,s_x) == IEOR(s_v0,s_v1) .and. &
IEOR(s_i,s_j) == IEOR(s_v1,s_v0)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W3(i,j,v1,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_v0) =  &
  W3_(s_v0, s_i)%array(i_v0, i_i)
end do
end do
! Z2 <-- T2(w,x,v0,v1) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_w, i_x) =  &
  T2_(s_v0, s_x, s_w)%array(i_v0, i_x, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_V, s_v0),&
                     - 4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,j)  <-- Z3
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
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccoo_ccvv_no3_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccoo_ccvv_no4_x0_type1_eri_v &
  (sj, ij, sv1, iv1, V2, W4, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv1, iv1
real(kind=8), intent(inout) :: V2(*), W4(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sj,sv1)

call set_symblock_Xav(sleft, W4, nir, nsym, psym) ! -> Xav (allocate) 
call g_sigma_ccoo_ccvv_no4_x0_type1_eri_v &
  (sj, ij, sv1, iv1, h2_i, Xav, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xav)

end subroutine g_if_sigma_ccoo_ccvv_no4_x0_type1_eri_v



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
subroutine g_sigma_ccoo_ccvv_no4_x0_type1_eri_v &
  (s_j, i_j, s_v1, i_v1, V2_, W4_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W4_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_a0, i_a0, s_i, i_i
! W4(i,j,v1,v0) += (    1.00000000) V2(v1,j,v0,a0) D1(i,a0) 
do s_v0 = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_i,s_j) == IEOR(s_v1,s_v0) .and. & 
IEOR(s_v1,s_j) == IEOR(s_v0,s_a0) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- V2(v1,j,v0,a0) 
allocate(Z1_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z1_(i_v0, i_a0) =  &
  V2_(s_a0, s_v0, s_j)%array(i_a0, i_v0, i_j)
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

! Z3 <-- W4(i,j,v1,v0) 
allocate(Z3_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_v0),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_v0))

! W4(i,j,v1,v0)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
W4_(s_v0, s_i)%array(i_v0, i_i) = &
    W4_(s_v0, s_i)%array(i_v0, i_i) &
  + Z3_(i_v0, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_v0) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccoo_ccvv_no4_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccoo_ccvv_no4_x1_type1_eri_v &
  (sj, ij, sv1, iv1, T2, W4, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv1, iv1
real(kind=8), intent(inout) :: T2(*), W4(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sv1)

call set_symblock_Xav(sleft, W4, nir, nsym, psym) ! -> Xav (allocate) 
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccoo_ccvv_no4_x1_type1_eri_v &
  (sj, ij, sv1, iv1, av2_i, Xav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xav)

end subroutine g_if_sigma_ccoo_ccvv_no4_x1_type1_eri_v



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
subroutine g_sigma_ccoo_ccvv_no4_x1_type1_eri_v &
  (s_j, i_j, s_v1, i_v1, T2_, W4_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W4_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_x, i_x, s_w, i_w, s_v0, i_v0, s_i, i_i
! S2(w,x,i,j) += (    2.00000000) T2(x,w,v0,v1) W4(i,j,v1,v0) 
do s_x = 0, nir-1
do s_w = 0, nir-1
do s_v0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_j) .and. & 
IEOR(s_x,s_w) == IEOR(s_v0,s_v1) .and. &
IEOR(s_i,s_j) == IEOR(s_v1,s_v0)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W4(i,j,v1,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_v0) =  &
  W4_(s_v0, s_i)%array(i_v0, i_i)
end do
end do
! Z2 <-- T2(x,w,v0,v1) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_x, i_w) =  &
  T2_(s_v0, s_w, s_x)%array(i_v0, i_w, i_x)
end do
end do
end do

! Z3 <-- S2(w,x,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_x)*psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_V, s_v0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,j)  <-- Z3
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
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccoo_ccvv_no4_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccoo_ccvv_no5_x0_type1_eri_v &
  (sj, ij, sv1, iv1, T2, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccoo_ccvv_no5_x0_type1_eri_v &
  (sj, ij, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_ccoo_ccvv_no5_x0_type1_eri_v



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
subroutine g_sigma_ccoo_ccvv_no5_x0_type1_eri_v &
  (s_j, i_j, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_j, s_j
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

integer :: s_v0, i_v0, s_i, i_i, s_w, i_w, s_x, i_x
! S2(w,x,i,j) += (    8.00000000) V2(v1,j,v0,i) T2(w,x,v0,v1) 
do s_v0 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_j) .and. & 
IEOR(s_v1,s_j) == IEOR(s_v0,s_i) .and. &
IEOR(s_w,s_x) == IEOR(s_v0,s_v1)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- V2(v1,j,v0,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_v0) =  &
  V2_(s_i, s_v0, s_j)%array(i_i, i_v0, i_j)
end do
end do
! Z2 <-- T2(w,x,v0,v1) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_w, i_x) =  &
  T2_(s_v0, s_x, s_w)%array(i_v0, i_x, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_V, s_v0),&
                     8.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,j)  <-- Z3
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
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccoo_ccvv_no5_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ccoo_ccvv_no6_x0_type1_eri_v &
  (sj, ij, sv1, iv1, T2, V2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sj, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ccoo_ccvv_no6_x0_type1_eri_v &
  (sj, ij, sv1, iv1, av2_i, h2_i, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_ccoo_ccvv_no6_x0_type1_eri_v



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
subroutine g_sigma_ccoo_ccvv_no6_x0_type1_eri_v &
  (s_j, i_j, s_v1, i_v1, T2_, V2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_j, s_j
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

integer :: s_i, i_i, s_v0, i_v0, s_w, i_w, s_x, i_x
! S2(w,x,i,j) += (   -4.00000000) V2(v1,i,v0,j) T2(w,x,v0,v1) 
do s_i = 0, nir-1
do s_v0 = 0, nir-1
do s_w = 0, nir-1
do s_x = 0, nir-1
if( &
IEOR(s_w,s_x) == IEOR(s_i,s_j) .and. & 
IEOR(s_v1,s_i) == IEOR(s_v0,s_j) .and. &
IEOR(s_w,s_x) == IEOR(s_v0,s_v1)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- V2(v1,i,v0,j) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_v0) =  &
  V2_(s_j, s_v0, s_i)%array(i_j, i_v0, i_i)
end do
end do
! Z2 <-- T2(w,x,v0,v1) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))
do i_x = psym(I_BEGIN, I_C, s_x), psym(I_END, I_C, s_x)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_w, i_x) =  &
  T2_(s_v0, s_x, s_w)%array(i_v0, i_x, i_w)
end do
end do
end do

! Z3 <-- S2(w,x,i,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_C, s_x):psym(I_END,I_C, s_x)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_C, s_x),&
                     psym(I_LENGTH,I_V, s_v0),&
                     - 4.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,x,i,j)  <-- Z3
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
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_ccoo_ccvv_no6_x0_type1_eri_v

