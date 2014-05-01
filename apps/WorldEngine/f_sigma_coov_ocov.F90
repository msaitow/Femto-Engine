#include <sci/icmr/fsrc/f_mr.fh>


!       #                #########      #     #   # 
!  ########## ##########         #   #######  #   # 
!      #    #         #          #    # #     #   # 
!      #    #        #   ########     # #     #   # 
!     #     #     # #           #  ##########    #  
!    #   # #       #            #       #       #   
!   #     #         #    ########       #     ##    

!                                    Generated date : Sun Apr 20 10:26:18 2014



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no0_x0_type0_noeri &
  (sa, ia, Fc0, T2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Fc0
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no0_x0_type0_noeri &
  (sa, ia, Fc0, av2_i, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_coov_ocov_no0_x0_type0_noeri



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
subroutine g_sigma_coov_ocov_no0_x0_type0_noeri &
  (s_a, i_a, Fc0, T2_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Fc0
! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_w, i_w, s_a0, i_a0, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) Fc0 T2(a1,w,a0,a) D2(j,i,a1,a0) 
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a1,s_w) == IEOR(s_a0,s_a) .and. &
IEOR(s_j,s_i) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(j,i,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a1, i_a0) =  &
  D2_(s_a0, s_a1, s_i, s_j)%array(i_a0, i_a1, i_i, i_j)
end do
end do
end do
end do
! Z2 <-- T2(a1,w,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_w) =  &
  T2_(s_a0, s_w, s_a1)%array(i_a0, i_w, i_a1)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00*Fc0, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_coov_ocov_no0_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no1_x0_type0_noeri &
  (sa, ia, Fc0, T2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Fc0
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no1_x0_type0_noeri &
  (sa, ia, Fc0, av2_i, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_coov_ocov_no1_x0_type0_noeri



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
subroutine g_sigma_coov_ocov_no1_x0_type0_noeri &
  (s_a, i_a, Fc0, T2_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
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

integer :: s_a0, i_a0, s_w, i_w, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) Fc0 T2(a0,w,j,a) D1(i,a0) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a0,s_w) == IEOR(s_j,s_a) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,w,j,a) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_j, i_a0) =  &
  T2_(s_j, s_w, s_a0)%array(i_j, i_w, i_a0)
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

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00*Fc0, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_w, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no1_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no2_x0_type0_noeri &
  (sa, ia, T2, W0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), W0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcaa(sleft, W0, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_coov_ocov_no2_x0_type0_noeri &
  (sa, ia, av2_i, Xcaa, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no2_x0_type0_noeri



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
subroutine g_sigma_coov_ocov_no2_x0_type0_noeri &
  (s_a, i_a, T2_, W0_, D2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_c0, i_c0, s_a0, i_a0, s_j, i_j, s_i, i_i
! W0(c0,j,i,a) += (    1.00000000) T2(a1,c0,a0,a) D2(j,i,a1,a0) 
do s_a1 = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_c0,s_j) == IEOR(s_i,s_a) .and. & 
IEOR(s_a1,s_c0) == IEOR(s_a0,s_a) .and. &
IEOR(s_j,s_i) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(j,i,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a1, i_a0) =  &
  D2_(s_a0, s_a1, s_i, s_j)%array(i_a0, i_a1, i_i, i_j)
end do
end do
end do
end do
! Z2 <-- T2(a1,c0,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_c0) =  &
  T2_(s_a0, s_c0, s_a1)%array(i_a0, i_c0, i_a1)
end do
end do
end do

! Z3 <-- W0(c0,j,i,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! W0(c0,j,i,a)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W0_(s_i, s_j, s_c0)%array(i_i, i_j, i_c0) = &
    W0_(s_i, s_j, s_c0)%array(i_i, i_j, i_c0) &
  + Z3_(i_j, i_i, i_c0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_coov_ocov_no2_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no2_x1_type0_noeri &
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

sleft = sa

call set_symblock_Xcaa(sleft, W0, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no2_x1_type0_noeri &
  (sa, ia, Xcaa, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no2_x1_type0_noeri



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
subroutine g_sigma_coov_ocov_no2_x1_type0_noeri &
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
type(symblock3), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (    1.00000000) Fc1(w,c0) W0(c0,j,i,a) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_w,s_c0) == 0 .and. &
IEOR(s_c0,s_j) == IEOR(s_i,s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- W0(c0,j,i,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_c0) =  &
  W0_(s_i, s_j, s_c0)%array(i_i, i_j, i_c0)
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

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_coov_ocov_no2_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no3_x0_type0_noeri &
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
call set_symblock_Xaaaa(sleft, W1, nir, nsym, psym) ! -> Xaaaa (allocate) 
call g_sigma_coov_ocov_no3_x0_type0_noeri &
  (Xaaaa, d3, fc1, nir, nsym, psym, flops)

deallocate(Xaaaa)

end subroutine g_if_sigma_coov_ocov_no3_x0_type0_noeri



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
! RDM is rotated :: D3(j,i,a3,a2,a1,a0)  >> D3(a2,a3,a0,a1,i,j) 
! rowInd : @[j, "active"] @[i, "active"] @[a1, "active"] @[a0, "active"] 
! summedInd : @[a3, "active"] @[a2, "active"] 
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_coov_ocov_no3_x0_type0_noeri &
  (W1_, D3_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W1_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_i, i_i, s_a3, i_a3, s_a2, i_a2, s_a1, i_a1
integer :: s_a0, i_a0
! W1(j,a1,i,a0) += (    1.00000000) D3(j,i,a3,a2,a1,a0) Fc1(a3,a2) 
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a3 = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_j,s_a1) == IEOR(s_i,s_a0) .and. & 
IEOR(IEOR(s_j,s_i),s_a3) == IEOR(IEOR(s_a2,s_a1),s_a0) .and. &
IEOR(s_a3,s_a2) == 0) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z2 <-- Fc1(a3,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a2) =  &
  Fc1_(s_a2, s_a3)%array(i_a2, i_a3)
end do
end do

! Z3 <-- W1(j,a1,i,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     D3_(s_j, s_i, s_a1, s_a0, s_a3, s_a2)%array,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0))

! W1(j,a1,i,a0)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W1_(s_a0, s_i, s_a1, s_j)%array(i_a0, i_i, i_a1, i_j) = &
    W1_(s_a0, s_i, s_a1, s_j)%array(i_a0, i_i, i_a1, i_j) &
  + Z3_(i_j, i_i, i_a1, i_a0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

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

end subroutine g_sigma_coov_ocov_no3_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no3_x1_type0_noeri &
  (sa, ia, T2, W1, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), W1(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xaaaa(sleft, W1, nir, nsym, psym) ! -> Xaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no3_x1_type0_noeri &
  (sa, ia, av2_i, Xaaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaaaa)

end subroutine g_if_sigma_coov_ocov_no3_x1_type0_noeri



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
subroutine g_sigma_coov_ocov_no3_x1_type0_noeri &
  (s_a, i_a, T2_, W1_, S2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W1_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_w, i_w, s_a0, i_a0, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) T2(a1,w,a0,a) W1(j,a1,i,a0) 
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a1,s_w) == IEOR(s_a0,s_a) .and. &
IEOR(s_j,s_a1) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W1(j,a1,i,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a1, i_a0) =  &
  W1_(s_a0, s_i, s_a1, s_j)%array(i_a0, i_i, i_a1, i_j)
end do
end do
end do
end do
! Z2 <-- T2(a1,w,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_w) =  &
  T2_(s_a0, s_w, s_a1)%array(i_a0, i_w, i_a1)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_coov_ocov_no3_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no4_x0_type0_noeri &
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
call set_symblock_Xaaaa(sleft, W2, nir, nsym, psym) ! -> Xaaaa (allocate) 
call g_sigma_coov_ocov_no4_x0_type0_noeri &
  (Xaaaa, d2, fc1, nir, nsym, psym, flops)

deallocate(Xaaaa)

end subroutine g_if_sigma_coov_ocov_no4_x0_type0_noeri



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
subroutine g_sigma_coov_ocov_no4_x0_type0_noeri &
  (W2_, D2_, Fc1_, nir, nsym, psym, flops)

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
type(symblock4), intent(inout) :: W2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a2, i_a2, s_a1, i_a1, s_a0, i_a0, s_j, i_j
! W2(i,a1,a0,j) += (    1.00000000) D2(i,a2,a1,a0) Fc1(j,a2) 
do s_i = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_j = 0, nir-1
if( &
IEOR(s_i,s_a1) == IEOR(s_a0,s_j) .and. & 
IEOR(s_i,s_a2) == IEOR(s_a1,s_a0) .and. &
IEOR(s_j,s_a2) == 0) then

if(psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- Fc1(j,a2) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_a2) =  &
  Fc1_(s_a2, s_j)%array(i_a2, i_j)
end do
end do
! Z2 <-- D2(i,a2,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_i, i_a1, i_a0) =  &
  D2_(s_a0, s_a1, s_a2, s_i)%array(i_a0, i_a1, i_a2, i_i)
end do
end do
end do
end do

! Z3 <-- W2(i,a1,a0,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j))

! W2(i,a1,a0,j)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
W2_(s_j, s_a0, s_a1, s_i)%array(i_j, i_a0, i_a1, i_i) = &
    W2_(s_j, s_a0, s_a1, s_i)%array(i_j, i_a0, i_a1, i_i) &
  + Z3_(i_j, i_i, i_a1, i_a0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no4_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no4_x1_type0_noeri &
  (sa, ia, T2, W2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), W2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xaaaa(sleft, W2, nir, nsym, psym) ! -> Xaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no4_x1_type0_noeri &
  (sa, ia, av2_i, Xaaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaaaa)

end subroutine g_if_sigma_coov_ocov_no4_x1_type0_noeri



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
subroutine g_sigma_coov_ocov_no4_x1_type0_noeri &
  (s_a, i_a, T2_, W2_, S2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_a1, i_a1, s_i, i_i, s_j, i_j
! S2(w,i,j,a) += (   -1.00000000) T2(a0,w,a1,a) W2(i,a1,a0,j) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_j = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a1,s_a) .and. &
IEOR(s_i,s_a1) == IEOR(s_a0,s_j)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W2(i,a1,a0,j) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_j, i_a1, i_a0) =  &
  W2_(s_j, s_a0, s_a1, s_i)%array(i_j, i_a0, i_a1, i_i)
end do
end do
end do
end do
! Z2 <-- T2(a0,w,a1,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_w) =  &
  T2_(s_a1, s_w, s_a0)%array(i_a1, i_w, i_a0)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_i, i_j, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no4_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no5_x0_type0_noeri &
  (sa, ia, T2, W3, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), W3(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcaa(sleft, W3, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_coov_ocov_no5_x0_type0_noeri &
  (sa, ia, av2_i, Xcaa, d1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no5_x0_type0_noeri



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
subroutine g_sigma_coov_ocov_no5_x0_type0_noeri &
  (s_a, i_a, T2_, W3_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W3_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_c0, i_c0, s_j, i_j, s_i, i_i
! W3(c0,i,j,a) += (    1.00000000) T2(a0,c0,j,a) D1(i,a0) 
do s_a0 = 0, nir-1
do s_c0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_c0,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a0,s_c0) == IEOR(s_j,s_a) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,c0,j,a) 
allocate(Z1_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z1_(i_c0, i_j, i_a0) =  &
  T2_(s_j, s_c0, s_a0)%array(i_j, i_c0, i_a0)
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

! Z3 <-- W3(c0,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_j))

! W3(c0,i,j,a)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
W3_(s_j, s_i, s_c0)%array(i_j, i_i, i_c0) = &
    W3_(s_j, s_i, s_c0)%array(i_j, i_i, i_c0) &
  + Z3_(i_c0, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no5_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no5_x1_type0_noeri &
  (sa, ia, W3, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: W3(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa

call set_symblock_Xcaa(sleft, W3, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no5_x1_type0_noeri &
  (sa, ia, Xcaa, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no5_x1_type0_noeri



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
subroutine g_sigma_coov_ocov_no5_x1_type0_noeri &
  (s_a, i_a, W3_, S2_, Fc1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W3_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_c0, i_c0, s_i, i_i, s_j, i_j
! S2(w,i,j,a) += (    1.00000000) Fc1(w,c0) W3(c0,i,j,a) 
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_i = 0, nir-1
do s_j = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_w,s_c0) == 0 .and. &
IEOR(s_c0,s_i) == IEOR(s_j,s_a)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- W3(c0,i,j,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_j, i_c0) =  &
  W3_(s_j, s_i, s_c0)%array(i_j, i_i, i_c0)
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

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_i, i_j, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no5_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no6_x0_type0_noeri &
  (W4, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

real(kind=8), intent(inout) :: W4(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xaa(sleft, W4, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_coov_ocov_no6_x0_type0_noeri &
  (Xaa, d2, fc1, nir, nsym, psym, flops)

deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no6_x0_type0_noeri



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
subroutine g_sigma_coov_ocov_no6_x0_type0_noeri &
  (W4_, D2_, Fc1_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W4_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_a2, i_a2, s_a1, i_a1
! W4(i,a0) += (    1.00000000) D2(i,a0,a2,a1) Fc1(a2,a1) 
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

! Z3 <-- W4(i,a0) 
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

! W4(i,a0)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W4_(s_a0, s_i)%array(i_a0, i_i) = &
    W4_(s_a0, s_i)%array(i_a0, i_i) &
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

end subroutine g_sigma_coov_ocov_no6_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no6_x1_type0_noeri &
  (sa, ia, T2, W4, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), W4(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xaa(sleft, W4, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no6_x1_type0_noeri &
  (sa, ia, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no6_x1_type0_noeri



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
subroutine g_sigma_coov_ocov_no6_x1_type0_noeri &
  (s_a, i_a, T2_, W4_, S2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W4_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) T2(a0,w,j,a) W4(i,a0) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a0,s_w) == IEOR(s_j,s_a) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,w,j,a) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_j, i_a0) =  &
  T2_(s_j, s_w, s_a0)%array(i_j, i_w, i_a0)
end do
end do
end do
! Z2 <-- W4(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  W4_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_w, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no6_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no7_x0_type0_noeri &
  (W5, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

real(kind=8), intent(inout) :: W5(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xaaaa(sleft, W5, nir, nsym, psym) ! -> Xaaaa (allocate) 
call g_sigma_coov_ocov_no7_x0_type0_noeri &
  (Xaaaa, d2, fc1, nir, nsym, psym, flops)

deallocate(Xaaaa)

end subroutine g_if_sigma_coov_ocov_no7_x0_type0_noeri



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
subroutine g_sigma_coov_ocov_no7_x0_type0_noeri &
  (W5_, D2_, Fc1_, nir, nsym, psym, flops)

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
type(symblock4), intent(inout) :: W5_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_i, i_i, s_a0, i_a0, s_a2, i_a2, s_a1, i_a1
! W5(j,a0,i,a1) += (    1.00000000) D2(j,i,a0,a2) Fc1(a2,a1) 
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(s_j,s_a0) == IEOR(s_i,s_a1) .and. & 
IEOR(s_j,s_i) == IEOR(s_a0,s_a2) .and. &
IEOR(s_a2,s_a1) == 0) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- Fc1(a2,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_a2) =  &
  Fc1_(s_a1, s_a2)%array(i_a1, i_a2)
end do
end do
! Z2 <-- D2(j,i,a0,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_j, i_i, i_a0) =  &
  D2_(s_a2, s_a0, s_i, s_j)%array(i_a2, i_a0, i_i, i_j)
end do
end do
end do
end do

! Z3 <-- W5(j,a0,i,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W5(j,a0,i,a1)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W5_(s_a1, s_i, s_a0, s_j)%array(i_a1, i_i, i_a0, i_j) = &
    W5_(s_a1, s_i, s_a0, s_j)%array(i_a1, i_i, i_a0, i_j) &
  + Z3_(i_a1, i_j, i_i, i_a0)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no7_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no7_x1_type0_noeri &
  (sa, ia, T2, W5, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), W5(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xaaaa(sleft, W5, nir, nsym, psym) ! -> Xaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no7_x1_type0_noeri &
  (sa, ia, av2_i, Xaaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaaaa)

end subroutine g_if_sigma_coov_ocov_no7_x1_type0_noeri



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
subroutine g_sigma_coov_ocov_no7_x1_type0_noeri &
  (s_a, i_a, T2_, W5_, S2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W5_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_a1, i_a1, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) T2(a0,w,a1,a) W5(j,a0,i,a1) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a1,s_a) .and. &
IEOR(s_j,s_a0) == IEOR(s_i,s_a1)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- W5(j,a0,i,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a0, i_a1) =  &
  W5_(s_a1, s_i, s_a0, s_j)%array(i_a1, i_i, i_a0, i_j)
end do
end do
end do
end do
! Z2 <-- T2(a0,w,a1,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a1, i_w) =  &
  T2_(s_a1, s_w, s_a0)%array(i_a1, i_w, i_a0)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
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

end subroutine g_sigma_coov_ocov_no7_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no8_x0_type0_noeri &
  (sa, ia, T2, W6, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), W6(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcaa(sleft, W6, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_coov_ocov_no8_x0_type0_noeri &
  (sa, ia, av2_i, Xcaa, d1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no8_x0_type0_noeri



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
subroutine g_sigma_coov_ocov_no8_x0_type0_noeri &
  (s_a, i_a, T2_, W6_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W6_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_a1, i_a1, s_i, i_i
! W6(w,i,a1,a) += (    1.00000000) T2(a0,w,a1,a) D1(i,a0) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a1,s_a) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a1,s_a) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,w,a1,a) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_a1, i_a0) =  &
  T2_(s_a1, s_w, s_a0)%array(i_a1, i_w, i_a0)
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

! Z3 <-- W6(w,i,a1,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1))

! W6(w,i,a1,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W6_(s_a1, s_i, s_w)%array(i_a1, i_i, i_w) = &
    W6_(s_a1, s_i, s_w)%array(i_a1, i_i, i_w) &
  + Z3_(i_w, i_a1, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1) * &
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

end subroutine g_sigma_coov_ocov_no8_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no8_x1_type0_noeri &
  (sa, ia, W6, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: W6(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa

call set_symblock_Xcaa(sleft, W6, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no8_x1_type0_noeri &
  (sa, ia, Xcaa, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no8_x1_type0_noeri



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
subroutine g_sigma_coov_ocov_no8_x1_type0_noeri &
  (s_a, i_a, W6_, S2_, Fc1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W6_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_a1, i_a1, s_w, i_w, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) Fc1(j,a1) W6(w,i,a1,a) 
do s_j = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_j,s_a1) == 0 .and. &
IEOR(s_w,s_i) == IEOR(s_a1,s_a)) then

if(psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- Fc1(j,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_a1) =  &
  Fc1_(s_a1, s_j)%array(i_a1, i_j)
end do
end do
! Z2 <-- W6(w,i,a1,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w, i_i) =  &
  W6_(s_a1, s_i, s_w)%array(i_a1, i_i, i_w)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a1),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_w, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no8_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no9_x0_type0_noeri &
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
call set_symblock_Xcaav(sleft, W7, nir, nsym, psym) ! -> Xcaav (allocate) 
call g_sigma_coov_ocov_no9_x0_type0_noeri &
  (sa0, ia0, av2_i, Xcaav, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcaav)

end subroutine g_if_sigma_coov_ocov_no9_x0_type0_noeri



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
subroutine g_sigma_coov_ocov_no9_x0_type0_noeri &
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

integer :: s_w, i_w, s_a1, i_a1, s_v0, i_v0, s_j, i_j, s_i, i_i
! W7(w,j,i,v0) += (    1.00000000) T2(w,a1,v0,a0) D2(j,i,a1,a0) 
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_v0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_j) == IEOR(s_i,s_v0) .and. & 
IEOR(s_w,s_a1) == IEOR(s_v0,s_a0) .and. &
IEOR(s_j,s_i) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- T2(w,a1,v0,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_v0, i_a1) =  &
  T2_(s_v0, s_a1, s_w)%array(i_v0, i_a1, i_w)
end do
end do
end do
! Z2 <-- D2(j,i,a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_j, i_i) =  &
  D2_(s_a0, s_a1, s_i, s_j)%array(i_a0, i_a1, i_i, i_j)
end do
end do
end do

! Z3 <-- W7(w,j,i,v0) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0))

! W7(w,j,i,v0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
W7_(s_v0, s_i, s_j, s_w)%array(i_v0, i_i, i_j, i_w) = &
    W7_(s_v0, s_i, s_j, s_w)%array(i_v0, i_i, i_j, i_w) &
  + Z3_(i_w, i_v0, i_j, i_i)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0) * &
                psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_coov_ocov_no9_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no9_x1_type0_noeri &
  (sa, ia, W7, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: W7(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_Xcaav(sleft, W7, nir, nsym, psym) ! -> Xcaav (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no9_x1_type0_noeri &
  (sa, ia, Xcaav, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaav)

end subroutine g_if_sigma_coov_ocov_no9_x1_type0_noeri



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
subroutine g_sigma_coov_ocov_no9_x1_type0_noeri &
  (s_a, i_a, W7_, S2_, Fc1_, nir, nsym, psym, flops)

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
type(symblock4), intent(inout) :: W7_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_w, i_w, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) Fc1(v0,a) W7(w,j,i,v0) 
do s_v0 = 0, nir-1
do s_w = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_v0,s_a) == 0 .and. &
IEOR(s_w,s_j) == IEOR(s_i,s_v0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W7(w,j,i,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_j, i_i, i_v0) =  &
  W7_(s_v0, s_i, s_j, s_w)%array(i_v0, i_i, i_j, i_w)
end do
end do
end do
end do
! Z2 <-- Fc1(v0,a) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0) =  &
  Fc1_(s_a, s_v0)%array(i_a, i_v0)
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_V, s_v0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_w, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
                1 * &
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no9_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no10_x0_type0_noeri &
  (sj, ij, T2, W8, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: T2(*), W8(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sj, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sj

call set_symblock_Xcav(sleft, W8, nir, nsym, psym) ! -> Xcav (allocate) 
call g_sigma_coov_ocov_no10_x0_type0_noeri &
  (sj, ij, av2_i, Xcav, d1, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcav)

end subroutine g_if_sigma_coov_ocov_no10_x0_type0_noeri



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
subroutine g_sigma_coov_ocov_no10_x0_type0_noeri &
  (s_j, i_j, T2_, W8_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W8_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_v0, i_v0, s_i, i_i
! W8(w,i,j,v0) += (    1.00000000) T2(w,a0,v0,j) D1(i,a0) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_v0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_v0) .and. & 
IEOR(s_w,s_a0) == IEOR(s_v0,s_j) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(w,a0,v0,j) 
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

! Z3 <-- W8(w,i,j,v0) 
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

! W8(w,i,j,v0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
W8_(s_v0, s_i, s_w)%array(i_v0, i_i, i_w) = &
    W8_(s_v0, s_i, s_w)%array(i_v0, i_i, i_w) &
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

end subroutine g_sigma_coov_ocov_no10_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no10_x1_type0_noeri &
  (sa, ia, sj, ij, W8, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sj, ij
real(kind=8), intent(inout) :: W8(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sj

call set_symblock_Xcav(sleft, W8, nir, nsym, psym) ! -> Xcav (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no10_x1_type0_noeri &
  (sa, ia, sj, ij, Xcav, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcav)

end subroutine g_if_sigma_coov_ocov_no10_x1_type0_noeri



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
subroutine g_sigma_coov_ocov_no10_x1_type0_noeri &
  (s_a, i_a, s_j, i_j, W8_, S2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W8_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_w, i_w, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) Fc1(v0,a) W8(w,i,j,v0) 
do s_v0 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_v0,s_a) == 0 .and. &
IEOR(s_w,s_i) == IEOR(s_j,s_v0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W8(w,i,j,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_i, i_v0) =  &
  W8_(s_v0, s_i, s_w)%array(i_v0, i_i, i_w)
end do
end do
end do
! Z2 <-- Fc1(v0,a) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0) =  &
  Fc1_(s_a, s_v0)%array(i_a, i_v0)
end do

! Z3 <-- S2(w,i,j,a) 
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

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
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

end subroutine g_sigma_coov_ocov_no10_x1_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no11_x0_type0_noeri &
  (sa, ia, T2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no11_x0_type0_noeri &
  (sa, ia, av2_i, av2_i2, c4f, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_coov_ocov_no11_x0_type0_noeri



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
subroutine g_sigma_coov_ocov_no11_x0_type0_noeri &
  (s_a, i_a, T2_, S2_, C4_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: C4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_w, i_w, s_a0, i_a0, s_i, i_i, s_j, i_j
! S2(w,i,j,a) += (   -0.50000000) T2(a1,w,a0,a) C4(a0,a1,i,j) 
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_j = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a1,s_w) == IEOR(s_a0,s_a) .and. &
IEOR(s_a0,s_a1) == IEOR(s_i,s_j)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- C4(a0,a1,i,j) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_j, i_a0, i_a1) =  &
  C4_(s_j, s_i, s_a1, s_a0)%array(i_j, i_i, i_a1, i_a0)
end do
end do
end do
end do
! Z2 <-- T2(a1,w,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a1, i_w) =  &
  T2_(s_a0, s_w, s_a1)%array(i_a0, i_w, i_a1)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     - 0.50000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_i, i_j, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_j) * &
                psym(I_LENGTH,I_C, s_w) * &
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

end subroutine g_sigma_coov_ocov_no11_x0_type0_noeri



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no0_x0_type1_eri_c &
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

call set_symblock_Xcaaaa(sleft, W9, nir, nsym, psym) ! -> Xcaaaa (allocate) 
call g_sigma_coov_ocov_no0_x0_type1_eri_c &
  (sw, iw, h2_i, Xcaaaa, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaaaa)

end subroutine g_if_sigma_coov_ocov_no0_x0_type1_eri_c



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
! RDM is rotated :: D3(j,a3,a2,i,a1,a0)  >> D3(i,a2,a0,a1,a3,j) 
! rowInd : @[j, "active"] @[i, "active"] @[a1, "active"] @[a0, "active"] 
! summedInd : @[a2, "active"] @[a3, "active"] 
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_coov_ocov_no0_x0_type1_eri_c &
  (s_w, i_w, V2_, W9_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: W9_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_c0, i_c0, s_a2, i_a2, s_j, i_j, s_i, i_i
integer :: s_a1, i_a1, s_a0, i_a0
! W9(w,c0,j,a1,i,a0) += (    1.00000000) V2(w,a3,c0,a2) D3(j,a3,a2,i,a1,a0) 
do s_a3 = 0, nir-1
do s_c0 = 0, nir-1
do s_a2 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(IEOR(s_w,s_c0),s_j) == IEOR(IEOR(s_a1,s_i),s_a0) .and. & 
IEOR(s_w,s_a3) == IEOR(s_c0,s_a2) .and. &
IEOR(IEOR(s_j,s_a3),s_a2) == IEOR(IEOR(s_i,s_a1),s_a0)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3) > 0) then

! Z1 <-- D3(i,a2,a0,a1,a3,j) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a1, i_a0, i_a2, i_a3) =  &
  D3_(s_j, s_a3, s_a1, s_a0, s_a2, s_i)%array(i_j, i_a3, i_a1, i_a0, i_a2, i_i)
end do
end do
end do
end do
end do
end do
! Z2 <-- V2(w,a3,c0,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a3, i_c0) =  &
  V2_(s_a2, s_c0, s_a3)%array(i_a2, i_c0, i_a3)
end do
end do
end do

! Z3 <-- W9(w,c0,j,a1,i,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0))

! W9(w,c0,j,a1,i,a0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W9_(s_a0, s_i, s_a1, s_j, s_c0)%array(i_a0, i_i, i_a1, i_j, i_c0) = &
    W9_(s_a0, s_i, s_a1, s_j, s_c0)%array(i_a0, i_i, i_a1, i_j, i_c0) &
  + Z3_(i_j, i_i, i_a1, i_a0, i_c0)
end do
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * &
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

end subroutine g_sigma_coov_ocov_no0_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no0_x1_type1_eri_c &
  (sa, ia, sw, iw, T2, W9, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: T2(*), W9(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xcaaaa(sleft, W9, nir, nsym, psym) ! -> Xcaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no0_x1_type1_eri_c &
  (sa, ia, sw, iw, av2_i, Xcaaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaaaa)

end subroutine g_if_sigma_coov_ocov_no0_x1_type1_eri_c



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
subroutine g_sigma_coov_ocov_no0_x1_type1_eri_c &
  (s_a, i_a, s_w, i_w, T2_, W9_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: W9_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_c0, i_c0, s_a0, i_a0, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (    1.00000000) T2(a1,c0,a0,a) W9(w,c0,j,a1,i,a0) 
do s_a1 = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a1,s_c0) == IEOR(s_a0,s_a) .and. &
IEOR(IEOR(s_w,s_c0),s_j) == IEOR(IEOR(s_a1,s_i),s_a0)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W9(w,c0,j,a1,i,a0) 
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
  W9_(s_a0, s_i, s_a1, s_j, s_c0)%array(i_a0, i_i, i_a1, i_j, i_c0)
end do
end do
end do
end do
end do
! Z2 <-- T2(a1,c0,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a1, i_a0) =  &
  T2_(s_a0, s_c0, s_a1)%array(i_a0, i_c0, i_a1)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
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

end subroutine g_sigma_coov_ocov_no0_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no1_x0_type1_eri_c &
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

call set_symblock_Xcaaaa(sleft, W10, nir, nsym, psym) ! -> Xcaaaa (allocate) 
call g_sigma_coov_ocov_no1_x0_type1_eri_c &
  (sw, iw, h2_i, Xcaaaa, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaaaa)

end subroutine g_if_sigma_coov_ocov_no1_x0_type1_eri_c



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
! RDM is rotated :: D3(j,i,a3,a2,a1,a0)  >> D3(a2,a3,a0,a1,i,j) 
! rowInd : @[j, "active"] @[i, "active"] @[a1, "active"] @[a0, "active"] 
! summedInd : @[a3, "active"] @[a2, "active"] 
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_coov_ocov_no1_x0_type1_eri_c &
  (s_w, i_w, V2_, W10_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: W10_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a3, i_a3, s_a2, i_a2, s_j, i_j, s_i, i_i
integer :: s_a1, i_a1, s_a0, i_a0
! W10(w,c0,j,a1,i,a0) += (    1.00000000) V2(w,c0,a3,a2) D3(j,i,a3,a2,a1,a0) 
do s_c0 = 0, nir-1
do s_a3 = 0, nir-1
do s_a2 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(IEOR(s_w,s_c0),s_j) == IEOR(IEOR(s_a1,s_i),s_a0) .and. & 
IEOR(s_w,s_c0) == IEOR(s_a3,s_a2) .and. &
IEOR(IEOR(s_j,s_i),s_a3) == IEOR(IEOR(s_a2,s_a1),s_a0)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z2 <-- V2(w,c0,a3,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a2, i_c0) =  &
  V2_(s_a2, s_a3, s_c0)%array(i_a2, i_a3, i_c0)
end do
end do
end do

! Z3 <-- W10(w,c0,j,a1,i,a0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0),&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     D3_(s_j, s_i, s_a1, s_a0, s_a3, s_a2)%array,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0))

! W10(w,c0,j,a1,i,a0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W10_(s_a0, s_i, s_a1, s_j, s_c0)%array(i_a0, i_i, i_a1, i_j, i_c0) = &
    W10_(s_a0, s_i, s_a1, s_j, s_c0)%array(i_a0, i_i, i_a1, i_j, i_c0) &
  + Z3_(i_j, i_i, i_a1, i_a0, i_c0)
end do
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0) * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z2_, Z3_)
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

end subroutine g_sigma_coov_ocov_no1_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no1_x1_type1_eri_c &
  (sa, ia, sw, iw, T2, W10, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: T2(*), W10(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xcaaaa(sleft, W10, nir, nsym, psym) ! -> Xcaaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no1_x1_type1_eri_c &
  (sa, ia, sw, iw, av2_i, Xcaaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaaaa)

end subroutine g_if_sigma_coov_ocov_no1_x1_type1_eri_c



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
subroutine g_sigma_coov_ocov_no1_x1_type1_eri_c &
  (s_a, i_a, s_w, i_w, T2_, W10_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_w, s_w
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: W10_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_c0, i_c0, s_a0, i_a0, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (    1.00000000) T2(a1,c0,a0,a) W10(w,c0,j,a1,i,a0) 
do s_a1 = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a1,s_c0) == IEOR(s_a0,s_a) .and. &
IEOR(IEOR(s_w,s_c0),s_j) == IEOR(IEOR(s_a1,s_i),s_a0)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W10(w,c0,j,a1,i,a0) 
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
  W10_(s_a0, s_i, s_a1, s_j, s_c0)%array(i_a0, i_i, i_a1, i_j, i_c0)
end do
end do
end do
end do
end do
! Z2 <-- T2(a1,c0,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a1, i_a0) =  &
  T2_(s_a0, s_c0, s_a1)%array(i_a0, i_c0, i_a1)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
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

end subroutine g_sigma_coov_ocov_no1_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no2_x0_type1_eri_c &
  (sa, ia, sc0, ic0, T2, W11, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sc0, ic0
real(kind=8), intent(inout) :: T2(*), W11(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sc0,sa)

call set_symblock_Xaa(sleft, W11, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_coov_ocov_no2_x0_type1_eri_c &
  (sa, ia, sc0, ic0, av2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no2_x0_type1_eri_c



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
subroutine g_sigma_coov_ocov_no2_x0_type1_eri_c &
  (s_a, i_a, s_c0, i_c0, T2_, W11_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c0, s_c0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W11_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_a1, i_a1, s_i, i_i, s_a2, i_a2
! W11(c0,i,a2,a) += (    1.00000000) T2(a0,c0,a1,a) D2(i,a2,a1,a0) 
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a2 = 0, nir-1
if( &
IEOR(s_c0,s_i) == IEOR(s_a2,s_a) .and. & 
IEOR(s_a0,s_c0) == IEOR(s_a1,s_a) .and. &
IEOR(s_i,s_a2) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(i,a2,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a2, i_a1, i_a0) =  &
  D2_(s_a0, s_a1, s_a2, s_i)%array(i_a0, i_a1, i_a2, i_i)
end do
end do
end do
end do
! Z2 <-- T2(a0,c0,a1,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  T2_(s_a1, s_c0, s_a0)%array(i_a1, i_c0, i_a0)
end do
end do

! Z3 <-- W11(c0,i,a2,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2))

! W11(c0,i,a2,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W11_(s_a2, s_i)%array(i_a2, i_i) = &
    W11_(s_a2, s_i)%array(i_a2, i_i) &
  + Z3_(i_i, i_a2)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) * &
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

end subroutine g_sigma_coov_ocov_no2_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no2_x1_type1_eri_c &
  (sa, ia, sc0, ic0, V2, W11, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sc0, ic0
real(kind=8), intent(inout) :: V2(*), W11(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sc0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sc0,sa)

call set_symblock_Xaa(sleft, W11, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no2_x1_type1_eri_c &
  (sa, ia, sc0, ic0, h2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no2_x1_type1_eri_c



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
subroutine g_sigma_coov_ocov_no2_x1_type1_eri_c &
  (s_a, i_a, s_c0, i_c0, V2_, W11_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c0, s_c0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W11_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_w, i_w, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (   -2.00000000) V2(c0,a2,w,j) W11(c0,i,a2,a) 
do s_a2 = 0, nir-1
do s_w = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_c0,s_a2) == IEOR(s_w,s_j) .and. &
IEOR(s_c0,s_i) == IEOR(s_a2,s_a)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(c0,a2,w,j) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_j, i_a2) =  &
  V2_(s_j, s_w, s_a2)%array(i_j, i_w, i_a2)
end do
end do
end do
! Z2 <-- W11(c0,i,a2,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_i) =  &
  W11_(s_a2, s_i)%array(i_a2, i_i)
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a2),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_w, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no2_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no3_x0_type1_eri_c &
  (sa, ia, sc0, ic0, T2, W12, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sc0, ic0
real(kind=8), intent(inout) :: T2(*), W12(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sc0,sa)

call set_symblock_Xaa(sleft, W12, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_coov_ocov_no3_x0_type1_eri_c &
  (sa, ia, sc0, ic0, av2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no3_x0_type1_eri_c



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
subroutine g_sigma_coov_ocov_no3_x0_type1_eri_c &
  (s_a, i_a, s_c0, i_c0, T2_, W12_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c0, s_c0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W12_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_a1, i_a1, s_i, i_i, s_a2, i_a2
! W12(c0,i,a2,a) += (    1.00000000) T2(a0,c0,a1,a) D2(i,a2,a1,a0) 
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a2 = 0, nir-1
if( &
IEOR(s_c0,s_i) == IEOR(s_a2,s_a) .and. & 
IEOR(s_a0,s_c0) == IEOR(s_a1,s_a) .and. &
IEOR(s_i,s_a2) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(i,a2,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a2, i_a1, i_a0) =  &
  D2_(s_a0, s_a1, s_a2, s_i)%array(i_a0, i_a1, i_a2, i_i)
end do
end do
end do
end do
! Z2 <-- T2(a0,c0,a1,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  T2_(s_a1, s_c0, s_a0)%array(i_a1, i_c0, i_a0)
end do
end do

! Z3 <-- W12(c0,i,a2,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2))

! W12(c0,i,a2,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W12_(s_a2, s_i)%array(i_a2, i_i) = &
    W12_(s_a2, s_i)%array(i_a2, i_i) &
  + Z3_(i_i, i_a2)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) * &
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

end subroutine g_sigma_coov_ocov_no3_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no3_x1_type1_eri_c &
  (sa, ia, sc0, ic0, V2, W12, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sc0, ic0
real(kind=8), intent(inout) :: V2(*), W12(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sc0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sc0,sa)

call set_symblock_Xaa(sleft, W12, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no3_x1_type1_eri_c &
  (sa, ia, sc0, ic0, h2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no3_x1_type1_eri_c



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
subroutine g_sigma_coov_ocov_no3_x1_type1_eri_c &
  (s_a, i_a, s_c0, i_c0, V2_, W12_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c0, s_c0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W12_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_j, i_j, s_a2, i_a2, s_i, i_i
! S2(w,i,j,a) += (    1.00000000) V2(c0,w,j,a2) W12(c0,i,a2,a) 
do s_w = 0, nir-1
do s_j = 0, nir-1
do s_a2 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_c0,s_w) == IEOR(s_j,s_a2) .and. &
IEOR(s_c0,s_i) == IEOR(s_a2,s_a)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(c0,w,j,a2) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_j, i_a2) =  &
  V2_(s_a2, s_j, s_w)%array(i_a2, i_j, i_w)
end do
end do
end do
! Z2 <-- W12(c0,i,a2,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_i) =  &
  W12_(s_a2, s_i)%array(i_a2, i_i)
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_w, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no3_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no4_x0_type1_eri_c &
  (sw, iw, V2, W14, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), W14(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_Xcaa(sleft, W14, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_coov_ocov_no4_x0_type1_eri_c &
  (sw, iw, h2_i, Xcaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no4_x0_type1_eri_c



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
subroutine g_sigma_coov_ocov_no4_x0_type1_eri_c &
  (s_w, i_w, V2_, W14_, D2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W14_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_c0, i_c0, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W14(w,c0,i,a0) += (    1.00000000) V2(w,a2,c0,a1) D2(i,a1,a2,a0) 
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

! Z3 <-- W14(w,c0,i,a0) 
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

! W14(w,c0,i,a0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W14_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0) = &
    W14_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0) &
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

end subroutine g_sigma_coov_ocov_no4_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no4_x1_type1_eri_c &
  (sa, ia, sw, iw, T2, W14, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: T2(*), W14(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xcaa(sleft, W14, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no4_x1_type1_eri_c &
  (sa, ia, sw, iw, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no4_x1_type1_eri_c



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
subroutine g_sigma_coov_ocov_no4_x1_type1_eri_c &
  (s_a, i_a, s_w, i_w, T2_, W14_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_w, s_w
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

integer :: s_a0, i_a0, s_c0, i_c0, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (    1.00000000) T2(a0,c0,j,a) W14(w,c0,i,a0) 
do s_a0 = 0, nir-1
do s_c0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a0,s_c0) == IEOR(s_j,s_a) .and. &
IEOR(s_w,s_c0) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(a0,c0,j,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_a0, i_c0) =  &
  T2_(s_j, s_c0, s_a0)%array(i_j, i_c0, i_a0)
end do
end do
end do
! Z2 <-- W14(w,c0,i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_c0, i_i) =  &
  W14_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no4_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no5_x0_type1_eri_c &
  (sw, iw, V2, W15, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), W15(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_Xcaa(sleft, W15, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_coov_ocov_no5_x0_type1_eri_c &
  (sw, iw, h2_i, Xcaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no5_x0_type1_eri_c



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
subroutine g_sigma_coov_ocov_no5_x0_type1_eri_c &
  (s_w, i_w, V2_, W15_, D2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W15_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_a2, i_a2, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W15(w,c0,i,a0) += (    1.00000000) V2(w,c0,a2,a1) D2(i,a0,a2,a1) 
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

! Z3 <-- W15(w,c0,i,a0) 
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

! W15(w,c0,i,a0)  <-- Z3
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W15_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0) = &
    W15_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0) &
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

end subroutine g_sigma_coov_ocov_no5_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no5_x1_type1_eri_c &
  (sa, ia, sw, iw, T2, W15, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: T2(*), W15(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xcaa(sleft, W15, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no5_x1_type1_eri_c &
  (sa, ia, sw, iw, av2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no5_x1_type1_eri_c



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
subroutine g_sigma_coov_ocov_no5_x1_type1_eri_c &
  (s_a, i_a, s_w, i_w, T2_, W15_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_w, s_w
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

integer :: s_a0, i_a0, s_c0, i_c0, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (    1.00000000) T2(a0,c0,j,a) W15(w,c0,i,a0) 
do s_a0 = 0, nir-1
do s_c0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a0,s_c0) == IEOR(s_j,s_a) .and. &
IEOR(s_w,s_c0) == IEOR(s_i,s_a0)) then

if(psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(a0,c0,j,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_a0, i_c0) =  &
  T2_(s_j, s_c0, s_a0)%array(i_j, i_c0, i_a0)
end do
end do
end do
! Z2 <-- W15(w,c0,i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_c0, i_i) =  &
  W15_(s_a0, s_i, s_c0)%array(i_a0, i_i, i_c0)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_C, s_c0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no5_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no6_x0_type1_eri_c &
  (sa, ia, sw, iw, T2, V2, W17, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: T2(*), V2(*), W17(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sw,sa)

call set_symblock_Xaa(sleft, W17, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_coov_ocov_no6_x0_type1_eri_c &
  (sa, ia, sw, iw, av2_i, h2_i, Xaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no6_x0_type1_eri_c



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
subroutine g_sigma_coov_ocov_no6_x0_type1_eri_c &
  (s_a, i_a, s_w, i_w, T2_, V2_, W17_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W17_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_c0, i_c0, s_a1, i_a1, s_a0, i_a0
! W17(w,a0,a2,a) += (    1.00000000) V2(w,a2,c0,a1) T2(a0,c0,a1,a) 
do s_a2 = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_a2,s_a) .and. & 
IEOR(s_w,s_a2) == IEOR(s_c0,s_a1) .and. &
IEOR(s_a0,s_c0) == IEOR(s_a1,s_a)) then

if(psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(w,a2,c0,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_c0, i_a1) =  &
  V2_(s_a1, s_c0, s_a2)%array(i_a1, i_c0, i_a2)
end do
end do
end do
! Z2 <-- T2(a0,c0,a1,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a1, i_a0) =  &
  T2_(s_a1, s_c0, s_a0)%array(i_a1, i_c0, i_a0)
end do
end do
end do

! Z3 <-- W17(w,a0,a2,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a2))

! W17(w,a0,a2,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W17_(s_a2, s_a0)%array(i_a2, i_a0) = &
    W17_(s_a2, s_a0)%array(i_a2, i_a0) &
  + Z3_(i_a2, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a2) * &
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

end subroutine g_sigma_coov_ocov_no6_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no6_x1_type1_eri_c &
  (sa, ia, sw, iw, W17, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: W17(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sw,sa)

call set_symblock_Xaa(sleft, W17, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no6_x1_type1_eri_c &
  (sa, ia, sw, iw, Xaa, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no6_x1_type1_eri_c



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
subroutine g_sigma_coov_ocov_no6_x1_type1_eri_c &
  (s_a, i_a, s_w, i_w, W17_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W17_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_a2, i_a2, s_a0, i_a0, s_i, i_i
! S2(w,i,j,a) += (    1.00000000) D2(j,a2,a0,i) W17(w,a0,a2,a) 
do s_j = 0, nir-1
do s_a2 = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_j,s_a2) == IEOR(s_a0,s_i) .and. &
IEOR(s_w,s_a0) == IEOR(s_a2,s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(j,a2,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a2, i_a0) =  &
  D2_(s_i, s_a0, s_a2, s_j)%array(i_i, i_a0, i_a2, i_j)
end do
end do
end do
end do
! Z2 <-- W17(w,a0,a2,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a0) =  &
  W17_(s_a2, s_a0)%array(i_a2, i_a0)
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
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

end subroutine g_sigma_coov_ocov_no6_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no7_x0_type1_eri_c &
  (sa, ia, sw, iw, T2, V2, W18, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: T2(*), V2(*), W18(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sw,sa)

call set_symblock_Xaa(sleft, W18, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_coov_ocov_no7_x0_type1_eri_c &
  (sa, ia, sw, iw, av2_i, h2_i, Xaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no7_x0_type1_eri_c



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
subroutine g_sigma_coov_ocov_no7_x0_type1_eri_c &
  (s_a, i_a, s_w, i_w, T2_, V2_, W18_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
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

integer :: s_c0, i_c0, s_a2, i_a2, s_a1, i_a1, s_a0, i_a0
! W18(w,a0,a2,a) += (    1.00000000) V2(w,c0,a2,a1) T2(a0,c0,a1,a) 
do s_c0 = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_a2,s_a) .and. & 
IEOR(s_w,s_c0) == IEOR(s_a2,s_a1) .and. &
IEOR(s_a0,s_c0) == IEOR(s_a1,s_a)) then

if(psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(w,c0,a2,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_c0, i_a1) =  &
  V2_(s_a1, s_a2, s_c0)%array(i_a1, i_a2, i_c0)
end do
end do
end do
! Z2 <-- T2(a0,c0,a1,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a1, i_a0) =  &
  T2_(s_a1, s_c0, s_a0)%array(i_a1, i_c0, i_a0)
end do
end do
end do

! Z3 <-- W18(w,a0,a2,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a2))

! W18(w,a0,a2,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W18_(s_a2, s_a0)%array(i_a2, i_a0) = &
    W18_(s_a2, s_a0)%array(i_a2, i_a0) &
  + Z3_(i_a2, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a2) * &
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

end subroutine g_sigma_coov_ocov_no7_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no7_x1_type1_eri_c &
  (sa, ia, sw, iw, W18, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: W18(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sw,sa)

call set_symblock_Xaa(sleft, W18, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no7_x1_type1_eri_c &
  (sa, ia, sw, iw, Xaa, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no7_x1_type1_eri_c



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
subroutine g_sigma_coov_ocov_no7_x1_type1_eri_c &
  (s_a, i_a, s_w, i_w, W18_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W18_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_i, i_i, s_a0, i_a0, s_a2, i_a2
! S2(w,i,j,a) += (    1.00000000) D2(j,i,a0,a2) W18(w,a0,a2,a) 
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_a2 = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_j,s_i) == IEOR(s_a0,s_a2) .and. &
IEOR(s_w,s_a0) == IEOR(s_a2,s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- D2(j,i,a0,a2) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a0, i_a2) =  &
  D2_(s_a2, s_a0, s_i, s_j)%array(i_a2, i_a0, i_i, i_j)
end do
end do
end do
end do
! Z2 <-- W18(w,a0,a2,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a2) =  &
  W18_(s_a2, s_a0)%array(i_a2, i_a0)
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
                1 * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no7_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no8_x0_type1_eri_c &
  (sa, ia, sw, iw, T2, V2, W20, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: T2(*), V2(*), W20(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sw,sa)

call set_symblock_Xaa(sleft, W20, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_coov_ocov_no8_x0_type1_eri_c &
  (sa, ia, sw, iw, av2_i, h2_i, Xaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no8_x0_type1_eri_c



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
subroutine g_sigma_coov_ocov_no8_x0_type1_eri_c &
  (s_a, i_a, s_w, i_w, T2_, V2_, W20_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
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

integer :: s_j, i_j, s_c0, i_c0, s_a1, i_a1, s_a0, i_a0
! W20(w,a0,j,a) += (    1.00000000) V2(w,j,c0,a1) T2(a0,c0,a1,a) 
do s_j = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_j,s_a) .and. & 
IEOR(s_w,s_j) == IEOR(s_c0,s_a1) .and. &
IEOR(s_a0,s_c0) == IEOR(s_a1,s_a)) then

if(psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(w,j,c0,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_c0, i_a1) =  &
  V2_(s_a1, s_c0, s_j)%array(i_a1, i_c0, i_j)
end do
end do
end do
! Z2 <-- T2(a0,c0,a1,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a1, i_a0) =  &
  T2_(s_a1, s_c0, s_a0)%array(i_a1, i_c0, i_a0)
end do
end do
end do

! Z3 <-- W20(w,a0,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j))

! W20(w,a0,j,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
W20_(s_j, s_a0)%array(i_j, i_a0) = &
    W20_(s_j, s_a0)%array(i_j, i_a0) &
  + Z3_(i_j, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no8_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no8_x1_type1_eri_c &
  (sa, ia, sw, iw, W20, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: W20(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sw,sa)

call set_symblock_Xaa(sleft, W20, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no8_x1_type1_eri_c &
  (sa, ia, sw, iw, Xaa, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no8_x1_type1_eri_c



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
subroutine g_sigma_coov_ocov_no8_x1_type1_eri_c &
  (s_a, i_a, s_w, i_w, W20_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
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

integer :: s_i, i_i, s_a0, i_a0, s_j, i_j
! S2(w,i,j,a) += (   -2.00000000) D1(i,a0) W20(w,a0,j,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_j = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_j,s_a)) then

if(psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W20(w,a0,j,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_a0) =  &
  W20_(s_j, s_a0)%array(i_j, i_a0)
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

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no8_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no9_x0_type1_eri_c &
  (sa, ia, sw, iw, T2, V2, W22, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: T2(*), V2(*), W22(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sw,sa)

call set_symblock_Xaa(sleft, W22, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_coov_ocov_no9_x0_type1_eri_c &
  (sa, ia, sw, iw, av2_i, h2_i, Xaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no9_x0_type1_eri_c



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
subroutine g_sigma_coov_ocov_no9_x0_type1_eri_c &
  (s_a, i_a, s_w, i_w, T2_, V2_, W22_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
integer, intent(in) :: i_a, s_a
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

integer :: s_c0, i_c0, s_j, i_j, s_a1, i_a1, s_a0, i_a0
! W22(w,a0,j,a) += (    1.00000000) V2(w,c0,j,a1) T2(a0,c0,a1,a) 
do s_c0 = 0, nir-1
do s_j = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_j,s_a) .and. & 
IEOR(s_w,s_c0) == IEOR(s_j,s_a1) .and. &
IEOR(s_a0,s_c0) == IEOR(s_a1,s_a)) then

if(psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(w,c0,j,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_c0, i_a1) =  &
  V2_(s_a1, s_j, s_c0)%array(i_a1, i_j, i_c0)
end do
end do
end do
! Z2 <-- T2(a0,c0,a1,a) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z2_(i_c0, i_a1, i_a0) =  &
  T2_(s_a1, s_c0, s_a0)%array(i_a1, i_c0, i_a0)
end do
end do
end do

! Z3 <-- W22(w,a0,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j))

! W22(w,a0,j,a)  <-- Z3
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
W22_(s_j, s_a0)%array(i_j, i_a0) = &
    W22_(s_j, s_a0)%array(i_j, i_a0) &
  + Z3_(i_j, i_a0)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no9_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no9_x1_type1_eri_c &
  (sa, ia, sw, iw, W22, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: W22(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sw,sa)

call set_symblock_Xaa(sleft, W22, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no9_x1_type1_eri_c &
  (sa, ia, sw, iw, Xaa, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no9_x1_type1_eri_c



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
subroutine g_sigma_coov_ocov_no9_x1_type1_eri_c &
  (s_a, i_a, s_w, i_w, W22_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_w, s_w
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
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_j, i_j
! S2(w,i,j,a) += (    1.00000000) D1(i,a0) W22(w,a0,j,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_j = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_j,s_a)) then

if(psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W22(w,a0,j,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_a0) =  &
  W22_(s_j, s_a0)%array(i_j, i_a0)
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

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j) * &
                psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no9_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no0_x0_type1_eri_o &
  (sj, ij, V2, W13, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: V2(*), W13(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sj, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sj

call set_symblock_Xaaa(sleft, W13, nir, nsym, psym) ! -> Xaaa (allocate) 
call g_sigma_coov_ocov_no0_x0_type1_eri_o &
  (sj, ij, h2_i, Xaaa, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_coov_ocov_no0_x0_type1_eri_o



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
! RDM is rotated :: D3(i,a3,a4,a2,a1,a0)  >> D3(a2,a4,a3,i,a0,a1) 
! rowInd : @[a1, "active"] @[a0, "active"] @[i, "active"] 
! summedInd : @[a3, "active"] @[a4, "active"] @[a2, "active"] 
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_coov_ocov_no0_x0_type1_eri_o &
  (s_j, i_j, V2_, W13_, D3_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W13_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_a4, i_a4, s_a2, i_a2, s_i, i_i, s_a1, i_a1
integer :: s_a0, i_a0
! W13(i,a1,a0,j) += (    1.00000000) V2(j,a3,a4,a2) D3(i,a3,a4,a2,a1,a0) 
do s_a3 = 0, nir-1
do s_a4 = 0, nir-1
do s_a2 = 0, nir-1
do s_i = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a1) == IEOR(s_a0,s_j) .and. & 
IEOR(s_j,s_a3) == IEOR(s_a4,s_a2) .and. &
IEOR(IEOR(s_i,s_a3),s_a4) == IEOR(IEOR(s_a2,s_a1),s_a0)) then

if(psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z2 <-- V2(j,a3,a4,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a4):psym(I_END,I_O, s_a4), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a4 = psym(I_BEGIN, I_O, s_a4), psym(I_END, I_O, s_a4)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a4, i_a2) =  &
  V2_(s_a2, s_a4, s_a3)%array(i_a2, i_a4, i_a3)
end do
end do
end do

! Z3 <-- W13(i,a1,a0,j) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     D3_(s_a1, s_a0, s_i, s_a3, s_a4, s_a2)%array,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_i))

! W13(i,a1,a0,j)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W13_(s_a0, s_a1, s_i)%array(i_a0, i_a1, i_i) = &
    W13_(s_a0, s_a1, s_i)%array(i_a0, i_a1, i_i) &
  + Z3_(i_a1, i_a0, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_i) * &
                1 * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

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

end subroutine g_sigma_coov_ocov_no0_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no0_x1_type1_eri_o &
  (sa, ia, sj, ij, T2, W13, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sj, ij
real(kind=8), intent(inout) :: T2(*), W13(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sj

call set_symblock_Xaaa(sleft, W13, nir, nsym, psym) ! -> Xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no0_x1_type1_eri_o &
  (sa, ia, sj, ij, av2_i, Xaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_coov_ocov_no0_x1_type1_eri_o



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
subroutine g_sigma_coov_ocov_no0_x1_type1_eri_o &
  (s_a, i_a, s_j, i_j, T2_, W13_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_j, s_j
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

integer :: s_a0, i_a0, s_w, i_w, s_a1, i_a1, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) T2(a0,w,a1,a) W13(i,a1,a0,j) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a1,s_a) .and. &
IEOR(s_i,s_a1) == IEOR(s_a0,s_j)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W13(i,a1,a0,j) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a1, i_a0) =  &
  W13_(s_a0, s_a1, s_i)%array(i_a0, i_a1, i_i)
end do
end do
end do
! Z2 <-- T2(a0,w,a1,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_w) =  &
  T2_(s_a1, s_w, s_a0)%array(i_a1, i_w, i_a0)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_i, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no0_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no1_x0_type1_eri_o &
  (sa4, ia4, V2, W16, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa4, ia4
real(kind=8), intent(inout) :: V2(*), W16(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa4, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_Xaa(sleft, W16, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_coov_ocov_no1_x0_type1_eri_o &
  (sa4, ia4, h2_i, Xaa, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no1_x0_type1_eri_o



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
subroutine g_sigma_coov_ocov_no1_x0_type1_eri_o &
  (s_a4, i_a4, V2_, W16_, D3_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W16_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a3, i_a3, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W16(i,a0) += (    1.00000000) V2(a4,a2,a3,a1) D3(i,a0,a4,a2,a3,a1) 
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

! Z3 <-- W16(i,a0) 
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

! W16(i,a0)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W16_(s_a0, s_i)%array(i_a0, i_i) = &
    W16_(s_a0, s_i)%array(i_a0, i_i) &
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

end subroutine g_sigma_coov_ocov_no1_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no2_x0_type1_eri_o &
  (sa1, ia1, V2, W19, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa1, ia1
real(kind=8), intent(inout) :: V2(*), W19(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa1

call set_symblock_Xaaa(sleft, W19, nir, nsym, psym) ! -> Xaaa (allocate) 
call g_sigma_coov_ocov_no2_x0_type1_eri_o &
  (sa1, ia1, h2_i, Xaaa, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_coov_ocov_no2_x0_type1_eri_o



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
! RDM is rotated :: D3(j,i,a4,a2,a0,a3)  >> D3(a2,a4,a3,a0,i,j) 
! rowInd : @[j, "active"] @[i, "active"] @[a0, "active"] 
! summedInd : @[a3, "active"] @[a4, "active"] @[a2, "active"] 
! -- Check2 is skipped 
!    -- is4RDM.second false
!    -- allRDM.second false
!    -- precedence    -1
subroutine g_sigma_coov_ocov_no2_x0_type1_eri_o &
  (s_a1, i_a1, V2_, W19_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a1, s_a1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W19_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_a4, i_a4, s_a2, i_a2, s_j, i_j, s_i, i_i
integer :: s_a0, i_a0
! W19(j,a0,i,a1) += (    1.00000000) V2(a1,a3,a4,a2) D3(j,i,a4,a2,a0,a3) 
do s_a3 = 0, nir-1
do s_a4 = 0, nir-1
do s_a2 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_j,s_a0) == IEOR(s_i,s_a1) .and. & 
IEOR(s_a1,s_a3) == IEOR(s_a4,s_a2) .and. &
IEOR(IEOR(s_j,s_i),s_a4) == IEOR(IEOR(s_a2,s_a0),s_a3)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z2 <-- V2(a1,a3,a4,a2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a4):psym(I_END,I_O, s_a4), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a4 = psym(I_BEGIN, I_O, s_a4), psym(I_END, I_O, s_a4)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a4, i_a2) =  &
  V2_(s_a2, s_a4, s_a3)%array(i_a2, i_a4, i_a3)
end do
end do
end do

! Z3 <-- W19(j,a0,i,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     1,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     D3_(s_j, s_i, s_a0, s_a3, s_a4, s_a2)%array,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0))

! W19(j,a0,i,a1)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W19_(s_i, s_a0, s_j)%array(i_i, i_a0, i_j) = &
    W19_(s_i, s_a0, s_j)%array(i_i, i_a0, i_j) &
  + Z3_(i_j, i_i, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                1 * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a4)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

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

end subroutine g_sigma_coov_ocov_no2_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no2_x1_type1_eri_o &
  (sa, ia, sa1, ia1, T2, W19, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), W19(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa1

call set_symblock_Xaaa(sleft, W19, nir, nsym, psym) ! -> Xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no2_x1_type1_eri_o &
  (sa, ia, sa1, ia1, av2_i, Xaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_coov_ocov_no2_x1_type1_eri_o



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
subroutine g_sigma_coov_ocov_no2_x1_type1_eri_o &
  (s_a, i_a, s_a1, i_a1, T2_, W19_, S2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W19_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_a0, i_a0, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) T2(w,a0,a,a1) W19(j,a0,i,a1) 
do s_w = 0, nir-1
do s_a0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_w,s_a0) == IEOR(s_a,s_a1) .and. &
IEOR(s_j,s_a0) == IEOR(s_i,s_a1)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W19(j,a0,i,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a0) =  &
  W19_(s_i, s_a0, s_j)%array(i_i, i_a0, i_j)
end do
end do
end do
! Z2 <-- T2(w,a0,a,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  T2_(s_a, s_a0, s_w)%array(i_a, i_a0, i_w)
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_coov_ocov_no2_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no3_x0_type1_eri_o &
  (sj, ij, V2, W21, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: V2(*), W21(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sj, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sj

call set_symblock_Xaaa(sleft, W21, nir, nsym, psym) ! -> Xaaa (allocate) 
call g_sigma_coov_ocov_no3_x0_type1_eri_o &
  (sj, ij, h2_i, Xaaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_coov_ocov_no3_x0_type1_eri_o



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
subroutine g_sigma_coov_ocov_no3_x0_type1_eri_o &
  (s_j, i_j, V2_, W21_, D2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W21_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a3, i_a3, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W21(i,a0,j,a1) += (    1.00000000) V2(j,a2,a3,a1) D2(i,a2,a3,a0) 
do s_a2 = 0, nir-1
do s_a3 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_j,s_a1) .and. & 
IEOR(s_j,s_a2) == IEOR(s_a3,s_a1) .and. &
IEOR(s_i,s_a2) == IEOR(s_a3,s_a0)) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3) > 0) then

! Z1 <-- V2(j,a2,a3,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3)))
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_a2, i_a3) =  &
  V2_(s_a1, s_a3, s_a2)%array(i_a1, i_a3, i_a2)
end do
end do
end do
! Z2 <-- D2(i,a2,a3,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a3, i_i, i_a0) =  &
  D2_(s_a0, s_a3, s_a2, s_i)%array(i_a0, i_a3, i_a2, i_i)
end do
end do
end do
end do

! Z3 <-- W21(i,a0,j,a1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W21(i,a0,j,a1)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W21_(s_a1, s_a0, s_i)%array(i_a1, i_a0, i_i) = &
    W21_(s_a1, s_a0, s_i)%array(i_a1, i_a0, i_i) &
  + Z3_(i_a1, i_i, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no3_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no3_x1_type1_eri_o &
  (sa, ia, sj, ij, T2, W21, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sj, ij
real(kind=8), intent(inout) :: T2(*), W21(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sj

call set_symblock_Xaaa(sleft, W21, nir, nsym, psym) ! -> Xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no3_x1_type1_eri_o &
  (sa, ia, sj, ij, av2_i, Xaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_coov_ocov_no3_x1_type1_eri_o



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
subroutine g_sigma_coov_ocov_no3_x1_type1_eri_o &
  (s_a, i_a, s_j, i_j, T2_, W21_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_j, s_j
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

integer :: s_a0, i_a0, s_w, i_w, s_a1, i_a1, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) T2(a0,w,a1,a) W21(i,a0,j,a1) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a1,s_a) .and. &
IEOR(s_i,s_a0) == IEOR(s_j,s_a1)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- W21(i,a0,j,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a1) =  &
  W21_(s_a1, s_a0, s_i)%array(i_a1, i_a0, i_i)
end do
end do
end do
! Z2 <-- T2(a0,w,a1,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a1, i_w) =  &
  T2_(s_a1, s_w, s_a0)%array(i_a1, i_w, i_a0)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
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

end subroutine g_sigma_coov_ocov_no3_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no4_x0_type1_eri_o &
  (sj, ij, V2, W23, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sj, ij
real(kind=8), intent(inout) :: V2(*), W23(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sj, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sj

call set_symblock_Xaaa(sleft, W23, nir, nsym, psym) ! -> Xaaa (allocate) 
call g_sigma_coov_ocov_no4_x0_type1_eri_o &
  (sj, ij, h2_i, Xaaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_coov_ocov_no4_x0_type1_eri_o



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
subroutine g_sigma_coov_ocov_no4_x0_type1_eri_o &
  (s_j, i_j, V2_, W23_, D2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W23_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_a3, i_a3, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W23(i,a0,j,a2) += (    1.00000000) V2(j,a2,a3,a1) D2(i,a0,a3,a1) 
do s_a2 = 0, nir-1
do s_a3 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_j,s_a2) .and. & 
IEOR(s_j,s_a2) == IEOR(s_a3,s_a1) .and. &
IEOR(s_i,s_a0) == IEOR(s_a3,s_a1)) then

if(psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(j,a2,a3,a1) 
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
! Z2 <-- D2(i,a0,a3,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a1, i_i, i_a0) =  &
  D2_(s_a1, s_a3, s_a0, s_i)%array(i_a1, i_a3, i_a0, i_i)
end do
end do
end do
end do

! Z3 <-- W23(i,a0,j,a2) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a2))

! W23(i,a0,j,a2)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W23_(s_a2, s_a0, s_i)%array(i_a2, i_a0, i_i) = &
    W23_(s_a2, s_a0, s_i)%array(i_a2, i_a0, i_i) &
  + Z3_(i_a2, i_i, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a2) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no4_x0_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no4_x1_type1_eri_o &
  (sa, ia, sj, ij, T2, W23, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sj, ij
real(kind=8), intent(inout) :: T2(*), W23(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sj

call set_symblock_Xaaa(sleft, W23, nir, nsym, psym) ! -> Xaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no4_x1_type1_eri_o &
  (sa, ia, sj, ij, av2_i, Xaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaaa)

end subroutine g_if_sigma_coov_ocov_no4_x1_type1_eri_o



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
subroutine g_sigma_coov_ocov_no4_x1_type1_eri_o &
  (s_a, i_a, s_j, i_j, T2_, W23_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_j, s_j
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
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_a2, i_a2, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) T2(a0,w,a2,a) W23(i,a0,j,a2) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a2 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a2,s_a) .and. &
IEOR(s_i,s_a0) == IEOR(s_j,s_a2)) then

if(psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- W23(i,a0,j,a2) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a2) =  &
  W23_(s_a2, s_a0, s_i)%array(i_a2, i_a0, i_i)
end do
end do
end do
! Z2 <-- T2(a0,w,a2,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a2, i_w) =  &
  T2_(s_a2, s_w, s_a0)%array(i_a2, i_w, i_a0)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a2),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_i, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no4_x1_type1_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no0_x0_type2_eri_o &
  (sa, ia, T2, W16, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: T2(*), W16(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xaa(sleft, W16, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no0_x0_type2_eri_o &
  (sa, ia, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no0_x0_type2_eri_o



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
subroutine g_sigma_coov_ocov_no0_x0_type2_eri_o &
  (s_a, i_a, T2_, W16_, S2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W16_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (   -0.50000000) T2(a0,w,j,a) W16(i,a0) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a0,s_w) == IEOR(s_j,s_a) .and. &
IEOR(s_i,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,w,j,a) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_j, i_a0) =  &
  T2_(s_j, s_w, s_a0)%array(i_j, i_w, i_a0)
end do
end do
end do
! Z2 <-- W16(i,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  W16_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 0.50000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_w, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no0_x0_type2_eri_o



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no0_x0_type1_eri_v &
  (sa, ia, sa0, ia0, T2, V2, W24, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), W24(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xca(sleft, W24, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_coov_ocov_no0_x0_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_coov_ocov_no0_x0_type1_eri_v



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
subroutine g_sigma_coov_ocov_no0_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, W24_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W24_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_v0, i_v0, s_c0, i_c0, s_a1, i_a1
! W24(w,a1,a0,a) += (    1.00000000) V2(a,w,v0,c0) T2(c0,a1,v0,a0) 
do s_w = 0, nir-1
do s_v0 = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(s_w,s_a1) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_v0,s_c0) .and. &
IEOR(s_c0,s_a1) == IEOR(s_v0,s_a0)) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- T2(c0,a1,v0,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_c0, i_v0) =  &
  T2_(s_v0, s_a1, s_c0)%array(i_v0, i_a1, i_c0)
end do
end do
end do
! Z2 <-- V2(a,w,v0,c0) 
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

! Z3 <-- W24(w,a1,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W24(w,a1,a0,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W24_(s_a1, s_w)%array(i_a1, i_w) = &
    W24_(s_a1, s_w)%array(i_a1, i_w) &
  + Z3_(i_a1, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
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

end subroutine g_sigma_coov_ocov_no0_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no0_x1_type1_eri_v &
  (sa, ia, sa0, ia0, W24, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: W24(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa0,sa)

call set_symblock_Xca(sleft, W24, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no0_x1_type1_eri_v &
  (sa, ia, sa0, ia0, Xca, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_coov_ocov_no0_x1_type1_eri_v



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
subroutine g_sigma_coov_ocov_no0_x1_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, W24_, S2_, D2_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W24_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_i, i_i, s_a1, i_a1, s_w, i_w
! S2(w,i,j,a) += (   -2.00000000) D2(j,i,a1,a0) W24(w,a1,a0,a) 
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_j,s_i) == IEOR(s_a1,s_a0) .and. &
IEOR(s_w,s_a1) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(j,i,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a1) =  &
  D2_(s_a0, s_a1, s_i, s_j)%array(i_a0, i_a1, i_i, i_j)
end do
end do
end do
! Z2 <-- W24(w,a1,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w) =  &
  W24_(s_a1, s_w)%array(i_a1, i_w)
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_coov_ocov_no0_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no1_x0_type1_eri_v &
  (sa, ia, sa0, ia0, T2, V2, W25, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), W25(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xcaaa(sleft, W25, nir, nsym, psym) ! -> Xcaaa (allocate) 
call g_sigma_coov_ocov_no1_x0_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, h2_i, Xcaaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaaa)

end subroutine g_if_sigma_coov_ocov_no1_x0_type1_eri_v



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
subroutine g_sigma_coov_ocov_no1_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, W25_, nir, nsym, psym, flops)

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
type(symblock4), intent(inout) :: W25_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_v0, i_v0, s_a2, i_a2, s_w, i_w, s_a1, i_a1
! W25(w,a1,a0,a3,a2,a) += (    1.00000000) V2(a,a3,v0,a2) T2(w,a1,v0,a0) 
do s_a3 = 0, nir-1
do s_v0 = 0, nir-1
do s_a2 = 0, nir-1
do s_w = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(IEOR(s_w,s_a1),s_a0) == IEOR(IEOR(s_a3,s_a2),s_a) .and. & 
IEOR(s_a,s_a3) == IEOR(s_v0,s_a2) .and. &
IEOR(s_w,s_a1) == IEOR(s_v0,s_a0)) then

if(psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- V2(a,a3,v0,a2) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z1_(i_a3, i_a2, i_v0) =  &
  V2_(s_a2, s_v0, s_a3)%array(i_a2, i_v0, i_a3)
end do
end do
end do
! Z2 <-- T2(w,a1,v0,a0) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_w, i_a1) =  &
  T2_(s_v0, s_a1, s_w)%array(i_v0, i_a1, i_w)
end do
end do
end do

! Z3 <-- W25(w,a1,a0,a3,a2,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2))

! W25(w,a1,a0,a3,a2,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W25_(s_a2, s_a3, s_a1, s_w)%array(i_a2, i_a3, i_a1, i_w) = &
    W25_(s_a2, s_a3, s_a1, s_w)%array(i_a2, i_a3, i_a1, i_w) &
  + Z3_(i_a3, i_a2, i_w, i_a1)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no1_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no1_x1_type1_eri_v &
  (sa, ia, sa0, ia0, W25, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: W25(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa0,sa)

call set_symblock_Xcaaa(sleft, W25, nir, nsym, psym) ! -> Xcaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no1_x1_type1_eri_v &
  (sa, ia, sa0, ia0, Xcaaa, av2_i2, d3, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaaa)

end subroutine g_if_sigma_coov_ocov_no1_x1_type1_eri_v



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
subroutine g_sigma_coov_ocov_no1_x1_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, W25_, S2_, D3_, nir, nsym, psym, flops)

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
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W25_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_a2, i_a2, s_a3, i_a3, s_i, i_i, s_a1, i_a1
integer :: s_w, i_w
! S2(w,i,j,a) += (   -1.00000000) D3(j,a2,a3,i,a1,a0) W25(w,a1,a0,a3,a2,a) 
do s_j = 0, nir-1
do s_a2 = 0, nir-1
do s_a3 = 0, nir-1
do s_i = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(IEOR(s_j,s_a2),s_a3) == IEOR(IEOR(s_i,s_a1),s_a0) .and. &
IEOR(IEOR(s_w,s_a1),s_a0) == IEOR(IEOR(s_a3,s_a2),s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D3(j,a2,a3,i,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a2, i_a3, i_a1) =  &
  D3_(s_a0, s_a1, s_i, s_a3, s_a2, s_j)%array(i_a0, i_a1, i_i, i_a3, i_a2, i_j)
end do
end do
end do
end do
end do
! Z2 <-- W25(w,a1,a0,a3,a2,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a3, i_a1, i_w) =  &
  W25_(s_a2, s_a3, s_a1, s_w)%array(i_a2, i_a3, i_a1, i_w)
end do
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

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

end subroutine g_sigma_coov_ocov_no1_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no2_x0_type1_eri_v &
  (sa, ia, sa0, ia0, T2, V2, W26, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: T2(*), V2(*), W26(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sa)

call set_symblock_Xca(sleft, W26, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_coov_ocov_no2_x0_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_coov_ocov_no2_x0_type1_eri_v



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
subroutine g_sigma_coov_ocov_no2_x0_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, T2_, V2_, W26_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W26_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_w, i_w, s_c0, i_c0, s_a1, i_a1
! W26(w,a1,a0,a) += (    1.00000000) V2(a,v0,w,c0) T2(c0,a1,v0,a0) 
do s_v0 = 0, nir-1
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(s_w,s_a1) == IEOR(s_a0,s_a) .and. & 
IEOR(s_a,s_v0) == IEOR(s_w,s_c0) .and. &
IEOR(s_c0,s_a1) == IEOR(s_v0,s_a0)) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- T2(c0,a1,v0,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_c0, i_v0) =  &
  T2_(s_v0, s_a1, s_c0)%array(i_v0, i_a1, i_c0)
end do
end do
end do
! Z2 <-- V2(a,v0,w,c0) 
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

! Z3 <-- W26(w,a1,a0,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W26(w,a1,a0,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W26_(s_a1, s_w)%array(i_a1, i_w) = &
    W26_(s_a1, s_w)%array(i_a1, i_w) &
  + Z3_(i_a1, i_w)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
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

end subroutine g_sigma_coov_ocov_no2_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no2_x1_type1_eri_v &
  (sa, ia, sa0, ia0, W26, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa0, ia0
real(kind=8), intent(inout) :: W26(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sa0,sa)

call set_symblock_Xca(sleft, W26, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no2_x1_type1_eri_v &
  (sa, ia, sa0, ia0, Xca, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_coov_ocov_no2_x1_type1_eri_v



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
subroutine g_sigma_coov_ocov_no2_x1_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, W26_, S2_, D2_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W26_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_i, i_i, s_a1, i_a1, s_w, i_w
! S2(w,i,j,a) += (    1.00000000) D2(j,i,a1,a0) W26(w,a1,a0,a) 
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_j,s_i) == IEOR(s_a1,s_a0) .and. &
IEOR(s_w,s_a1) == IEOR(s_a0,s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(j,i,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a1) =  &
  D2_(s_a0, s_a1, s_i, s_j)%array(i_a0, i_a1, i_i, i_j)
end do
end do
end do
! Z2 <-- W26(w,a1,a0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_w) =  &
  W26_(s_a1, s_w)%array(i_a1, i_w)
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_coov_ocov_no2_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no3_x0_type1_eri_v &
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

call set_symblock_Xcaaa(sleft, W27, nir, nsym, psym) ! -> Xcaaa (allocate) 
call g_sigma_coov_ocov_no3_x0_type1_eri_v &
  (sa, ia, sa0, ia0, av2_i, h2_i, Xcaaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaaa)

end subroutine g_if_sigma_coov_ocov_no3_x0_type1_eri_v



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
subroutine g_sigma_coov_ocov_no3_x0_type1_eri_v &
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
type(symblock4), intent(inout) :: W27_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_a3, i_a3, s_a2, i_a2, s_w, i_w, s_a1, i_a1
! W27(w,a1,a0,a3,a2,a) += (    1.00000000) V2(a,v0,a3,a2) T2(w,a1,v0,a0) 
do s_v0 = 0, nir-1
do s_a3 = 0, nir-1
do s_a2 = 0, nir-1
do s_w = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(IEOR(s_w,s_a1),s_a0) == IEOR(IEOR(s_a3,s_a2),s_a) .and. & 
IEOR(s_a,s_v0) == IEOR(s_a3,s_a2) .and. &
IEOR(s_w,s_a1) == IEOR(s_v0,s_a0)) then

if(psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- V2(a,v0,a3,a2) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z1_(i_a3, i_a2, i_v0) =  &
  V2_(s_a2, s_a3, s_v0)%array(i_a2, i_a3, i_v0)
end do
end do
end do
! Z2 <-- T2(w,a1,v0,a0) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z2_(i_v0, i_w, i_a1) =  &
  T2_(s_v0, s_a1, s_w)%array(i_v0, i_a1, i_w)
end do
end do
end do

! Z3 <-- W27(w,a1,a0,a3,a2,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2))

! W27(w,a1,a0,a3,a2,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W27_(s_a2, s_a3, s_a1, s_w)%array(i_a2, i_a3, i_a1, i_w) = &
    W27_(s_a2, s_a3, s_a1, s_w)%array(i_a2, i_a3, i_a1, i_w) &
  + Z3_(i_a3, i_a2, i_w, i_a1)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a1) * &
                psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no3_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no3_x1_type1_eri_v &
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

call set_symblock_Xcaaa(sleft, W27, nir, nsym, psym) ! -> Xcaaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no3_x1_type1_eri_v &
  (sa, ia, sa0, ia0, Xcaaa, av2_i2, d3, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaaa)

end subroutine g_if_sigma_coov_ocov_no3_x1_type1_eri_v



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
subroutine g_sigma_coov_ocov_no3_x1_type1_eri_v &
  (s_a, i_a, s_a0, i_a0, W27_, S2_, D3_, nir, nsym, psym, flops)

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
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W27_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_i, i_i, s_a3, i_a3, s_a2, i_a2, s_a1, i_a1
integer :: s_w, i_w
! S2(w,i,j,a) += (   -1.00000000) D3(j,i,a3,a2,a1,a0) W27(w,a1,a0,a3,a2,a) 
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a3 = 0, nir-1
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(IEOR(s_j,s_i),s_a3) == IEOR(IEOR(s_a2,s_a1),s_a0) .and. &
IEOR(IEOR(s_w,s_a1),s_a0) == IEOR(IEOR(s_a3,s_a2),s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D3(j,i,a3,a2,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a3, i_a2, i_a1) =  &
  D3_(s_a0, s_a1, s_a2, s_a3, s_i, s_j)%array(i_a0, i_a1, i_a2, i_a3, i_i, i_j)
end do
end do
end do
end do
end do
! Z2 <-- W27(w,a1,a0,a3,a2,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a2, i_a1, i_w) =  &
  W27_(s_a2, s_a3, s_a1, s_w)%array(i_a2, i_a3, i_a1, i_w)
end do
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

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

end subroutine g_sigma_coov_ocov_no3_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no4_x0_type1_eri_v &
  (sv0, iv0, T2, W28, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), W28(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W28, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_coov_ocov_no4_x0_type1_eri_v &
  (sv0, iv0, av2_i, Xcaa, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no4_x0_type1_eri_v



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
subroutine g_sigma_coov_ocov_no4_x0_type1_eri_v &
  (s_v0, i_v0, T2_, W28_, D2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W28_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_a1, i_a1, s_i, i_i, s_a2, i_a2
! W28(w,i,a2,v0) += (    1.00000000) T2(a0,w,a1,v0) D2(i,a2,a1,a0) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a2 = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a2,s_v0) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a1,s_v0) .and. &
IEOR(s_i,s_a2) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(i,a2,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a2, i_a1, i_a0) =  &
  D2_(s_a0, s_a1, s_a2, s_i)%array(i_a0, i_a1, i_a2, i_i)
end do
end do
end do
end do
! Z2 <-- T2(a0,w,a1,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_w) =  &
  T2_(s_a1, s_w, s_a0)%array(i_a1, i_w, i_a0)
end do
end do
end do

! Z3 <-- W28(w,i,a2,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2))

! W28(w,i,a2,v0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W28_(s_a2, s_i, s_w)%array(i_a2, i_i, i_w) = &
    W28_(s_a2, s_i, s_w)%array(i_a2, i_i, i_w) &
  + Z3_(i_i, i_a2, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) * &
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

end subroutine g_sigma_coov_ocov_no4_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no4_x1_type1_eri_v &
  (sa, ia, sv0, iv0, V2, W28, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W28(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W28, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no4_x1_type1_eri_v &
  (sa, ia, sv0, iv0, h2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no4_x1_type1_eri_v



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
subroutine g_sigma_coov_ocov_no4_x1_type1_eri_v &
  (s_a, i_a, s_v0, i_v0, V2_, W28_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W28_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_a2, i_a2, s_w, i_w, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) V2(v0,a,j,a2) W28(w,i,a2,v0) 
do s_j = 0, nir-1
do s_a2 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_v0,s_a) == IEOR(s_j,s_a2) .and. &
IEOR(s_w,s_i) == IEOR(s_a2,s_v0)) then

if(psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(v0,a,j,a2) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_a2) =  &
  V2_(s_a2, s_j, s_a)%array(i_a2, i_j, i_a)
end do
end do
! Z2 <-- W28(w,i,a2,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_w, i_i) =  &
  W28_(s_a2, s_i, s_w)%array(i_a2, i_i, i_w)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a2),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_w, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no4_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no5_x0_type1_eri_v &
  (sa, ia, sj, ij, T2, V2, W29, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sj, ij
real(kind=8), intent(inout) :: T2(*), V2(*), W29(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sj, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sa)

call set_symblock_Xca(sleft, W29, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_coov_ocov_no5_x0_type1_eri_v &
  (sa, ia, sj, ij, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_coov_ocov_no5_x0_type1_eri_v



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
subroutine g_sigma_coov_ocov_no5_x0_type1_eri_v &
  (s_a, i_a, s_j, i_j, T2_, V2_, W29_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W29_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_w, i_w, s_v0, i_v0, s_c0, i_c0, s_a0, i_a0
! W29(w,a0,j,a) += (    1.00000000) V2(a,w,v0,c0) T2(c0,a0,v0,j) 
do s_w = 0, nir-1
do s_v0 = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_j,s_a) .and. & 
IEOR(s_a,s_w) == IEOR(s_v0,s_c0) .and. &
IEOR(s_c0,s_a0) == IEOR(s_v0,s_j)) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- T2(c0,a0,v0,j) 
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
! Z2 <-- V2(a,w,v0,c0) 
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

! Z3 <-- W29(w,a0,j,a) 
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

! W29(w,a0,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W29_(s_a0, s_w)%array(i_a0, i_w) = &
    W29_(s_a0, s_w)%array(i_a0, i_w) &
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

end subroutine g_sigma_coov_ocov_no5_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no5_x1_type1_eri_v &
  (sa, ia, sj, ij, W29, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sj, ij
real(kind=8), intent(inout) :: W29(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sj,sa)

call set_symblock_Xca(sleft, W29, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no5_x1_type1_eri_v &
  (sa, ia, sj, ij, Xca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_coov_ocov_no5_x1_type1_eri_v



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
subroutine g_sigma_coov_ocov_no5_x1_type1_eri_v &
  (s_a, i_a, s_j, i_j, W29_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W29_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_w, i_w
! S2(w,i,j,a) += (   -2.00000000) D1(i,a0) W29(w,a0,j,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_j,s_a)) then

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
! Z2 <-- W29(w,a0,j,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  W29_(s_a0, s_w)%array(i_a0, i_w)
end do
end do

! Z3 <-- S2(w,i,j,a) 
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

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
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

end subroutine g_sigma_coov_ocov_no5_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no6_x0_type1_eri_v &
  (sa, ia, sv0, iv0, V2, W30, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W30(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sv0,sa)

call set_symblock_Xaa(sleft, W30, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_coov_ocov_no6_x0_type1_eri_v &
  (sa, ia, sv0, iv0, h2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no6_x0_type1_eri_v



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
subroutine g_sigma_coov_ocov_no6_x0_type1_eri_v &
  (s_a, i_a, s_v0, i_v0, V2_, W30_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W30_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a2, i_a2, s_i, i_i, s_a0, i_a0
! W30(i,a0,v0,a) += (    1.00000000) V2(v0,a1,a2,a) D2(i,a2,a1,a0) 
do s_a1 = 0, nir-1
do s_a2 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_v0,s_a) .and. & 
IEOR(s_v0,s_a1) == IEOR(s_a2,s_a) .and. &
IEOR(s_i,s_a2) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- D2(i,a2,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a0, i_a2, i_a1) =  &
  D2_(s_a0, s_a1, s_a2, s_i)%array(i_a0, i_a1, i_a2, i_i)
end do
end do
end do
end do
! Z2 <-- V2(v0,a1,a2,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1) =  &
  V2_(s_a, s_a2, s_a1)%array(i_a, i_a2, i_a1)
end do
end do

! Z3 <-- W30(i,a0,v0,a) 
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

! W30(i,a0,v0,a)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W30_(s_a0, s_i)%array(i_a0, i_i) = &
    W30_(s_a0, s_i)%array(i_a0, i_i) &
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

end subroutine g_sigma_coov_ocov_no6_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no6_x1_type1_eri_v &
  (sa, ia, sv0, iv0, T2, W30, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W30(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sv0,sa)

call set_symblock_Xaa(sleft, W30, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no6_x1_type1_eri_v &
  (sa, ia, sv0, iv0, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no6_x1_type1_eri_v



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
subroutine g_sigma_coov_ocov_no6_x1_type1_eri_v &
  (s_a, i_a, s_v0, i_v0, T2_, W30_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W30_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) T2(a0,w,j,v0) W30(i,a0,v0,a) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a0,s_w) == IEOR(s_j,s_v0) .and. &
IEOR(s_i,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,w,j,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_j, i_a0) =  &
  T2_(s_j, s_w, s_a0)%array(i_j, i_w, i_a0)
end do
end do
end do
! Z2 <-- W30(i,a0,v0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  W30_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_w, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no6_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no7_x0_type1_eri_v &
  (sa, ia, sj, ij, T2, V2, W31, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sj, ij
real(kind=8), intent(inout) :: T2(*), V2(*), W31(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sj, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sj,sa)

call set_symblock_Xca(sleft, W31, nir, nsym, psym) ! -> Xca (allocate) 
call g_sigma_coov_ocov_no7_x0_type1_eri_v &
  (sa, ia, sj, ij, av2_i, h2_i, Xca, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xca)

end subroutine g_if_sigma_coov_ocov_no7_x0_type1_eri_v



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
subroutine g_sigma_coov_ocov_no7_x0_type1_eri_v &
  (s_a, i_a, s_j, i_j, T2_, V2_, W31_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_j, s_j
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W31_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_w, i_w, s_c0, i_c0, s_a0, i_a0
! W31(w,a0,j,a) += (    1.00000000) V2(a,v0,w,c0) T2(c0,a0,v0,j) 
do s_v0 = 0, nir-1
do s_w = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_j,s_a) .and. & 
IEOR(s_a,s_v0) == IEOR(s_w,s_c0) .and. &
IEOR(s_c0,s_a0) == IEOR(s_v0,s_j)) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- T2(c0,a0,v0,j) 
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
! Z2 <-- V2(a,v0,w,c0) 
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

! Z3 <-- W31(w,a0,j,a) 
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

! W31(w,a0,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W31_(s_a0, s_w)%array(i_a0, i_w) = &
    W31_(s_a0, s_w)%array(i_a0, i_w) &
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

end subroutine g_sigma_coov_ocov_no7_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no7_x1_type1_eri_v &
  (sa, ia, sj, ij, W31, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sj, ij
real(kind=8), intent(inout) :: W31(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = IEOR(sj,sa)

call set_symblock_Xca(sleft, W31, nir, nsym, psym) ! -> Xca (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no7_x1_type1_eri_v &
  (sa, ia, sj, ij, Xca, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xca)

end subroutine g_if_sigma_coov_ocov_no7_x1_type1_eri_v



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
subroutine g_sigma_coov_ocov_no7_x1_type1_eri_v &
  (s_a, i_a, s_j, i_j, W31_, S2_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_j, s_j
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock2), intent(inout) :: D1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W31_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_w, i_w
! S2(w,i,j,a) += (    1.00000000) D1(i,a0) W31(w,a0,j,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_j,s_a)) then

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
! Z2 <-- W31(w,a0,j,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_w) =  &
  W31_(s_a0, s_w)%array(i_a0, i_w)
end do
end do

! Z3 <-- S2(w,i,j,a) 
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

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
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

end subroutine g_sigma_coov_ocov_no7_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no8_x0_type1_eri_v &
  (sa, ia, sv0, iv0, V2, W32, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W32(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sv0,sa)

call set_symblock_Xaa(sleft, W32, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_coov_ocov_no8_x0_type1_eri_v &
  (sa, ia, sv0, iv0, h2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no8_x0_type1_eri_v



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
subroutine g_sigma_coov_ocov_no8_x0_type1_eri_v &
  (s_a, i_a, s_v0, i_v0, V2_, W32_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
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

integer :: s_a2, i_a2, s_a1, i_a1, s_i, i_i, s_a0, i_a0
! W32(i,a0,v0,a) += (    1.00000000) V2(v0,a,a2,a1) D2(i,a0,a2,a1) 
do s_a2 = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_i,s_a0) == IEOR(s_v0,s_a) .and. & 
IEOR(s_v0,s_a) == IEOR(s_a2,s_a1) .and. &
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
! Z2 <-- V2(v0,a,a2,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_a1) =  &
  V2_(s_a1, s_a2, s_a)%array(i_a1, i_a2, i_a)
end do
end do

! Z3 <-- W32(i,a0,v0,a) 
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

! W32(i,a0,v0,a)  <-- Z3
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
                psym(I_LENGTH,I_O, s_a2)*psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no8_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no8_x1_type1_eri_v &
  (sa, ia, sv0, iv0, T2, W32, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W32(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sv0,sa)

call set_symblock_Xaa(sleft, W32, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no8_x1_type1_eri_v &
  (sa, ia, sv0, iv0, av2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_ocov_no8_x1_type1_eri_v



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
subroutine g_sigma_coov_ocov_no8_x1_type1_eri_v &
  (s_a, i_a, s_v0, i_v0, T2_, W32_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
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
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (   -1.00000000) T2(a0,w,j,v0) W32(i,a0,v0,a) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a0,s_w) == IEOR(s_j,s_v0) .and. &
IEOR(s_i,s_a0) == IEOR(s_v0,s_a)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(a0,w,j,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_j, i_a0) =  &
  T2_(s_j, s_w, s_a0)%array(i_j, i_w, i_a0)
end do
end do
end do
! Z2 <-- W32(i,a0,v0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_i) =  &
  W32_(s_a0, s_i)%array(i_a0, i_i)
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_w, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no8_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no9_x0_type1_eri_v &
  (sv0, iv0, T2, W33, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), W33(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W33, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_coov_ocov_no9_x0_type1_eri_v &
  (sv0, iv0, av2_i, Xcaa, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no9_x0_type1_eri_v



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
subroutine g_sigma_coov_ocov_no9_x0_type1_eri_v &
  (s_v0, i_v0, T2_, W33_, D2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W33_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_w, i_w, s_a1, i_a1, s_i, i_i, s_a2, i_a2
! W33(w,i,a2,v0) += (    1.00000000) T2(a0,w,a1,v0) D2(i,a2,a1,a0) 
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a2 = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_a2,s_v0) .and. & 
IEOR(s_a0,s_w) == IEOR(s_a1,s_v0) .and. &
IEOR(s_i,s_a2) == IEOR(s_a1,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(i,a2,a1,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a2, i_a1, i_a0) =  &
  D2_(s_a0, s_a1, s_a2, s_i)%array(i_a0, i_a1, i_a2, i_i)
end do
end do
end do
end do
! Z2 <-- T2(a0,w,a1,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_w) =  &
  T2_(s_a1, s_w, s_a0)%array(i_a1, i_w, i_a0)
end do
end do
end do

! Z3 <-- W33(w,i,a2,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2))

! W33(w,i,a2,v0)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W33_(s_a2, s_i, s_w)%array(i_a2, i_i, i_w) = &
    W33_(s_a2, s_i, s_w)%array(i_a2, i_i, i_w) &
  + Z3_(i_i, i_a2, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) * &
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

end subroutine g_sigma_coov_ocov_no9_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no9_x1_type1_eri_v &
  (sa, ia, sv0, iv0, V2, W33, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W33(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xcaa(sleft, W33, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no9_x1_type1_eri_v &
  (sa, ia, sv0, iv0, h2_i, Xcaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no9_x1_type1_eri_v



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
subroutine g_sigma_coov_ocov_no9_x1_type1_eri_v &
  (s_a, i_a, s_v0, i_v0, V2_, W33_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W33_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_a2, i_a2, s_w, i_w, s_i, i_i
! S2(w,i,j,a) += (    2.00000000) V2(v0,j,a2,a) W33(w,i,a2,v0) 
do s_j = 0, nir-1
do s_a2 = 0, nir-1
do s_w = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_v0,s_j) == IEOR(s_a2,s_a) .and. &
IEOR(s_w,s_i) == IEOR(s_a2,s_v0)) then

if(psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(v0,j,a2,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_a2) =  &
  V2_(s_a, s_a2, s_j)%array(i_a, i_a2, i_j)
end do
end do
! Z2 <-- W33(w,i,a2,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_w, i_i) =  &
  W33_(s_a2, s_i, s_w)%array(i_a2, i_i, i_w)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a2),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_w, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j) * &
                psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no9_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no10_x0_type1_eri_v &
  (sa, ia, sa2, ia2, T2, V2, W34, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa2, ia2
real(kind=8), intent(inout) :: T2(*), V2(*), W34(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa2, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcaa(sleft, W34, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_coov_ocov_no10_x0_type1_eri_v &
  (sa, ia, sa2, ia2, av2_i, h2_i, Xcaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no10_x0_type1_eri_v



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
subroutine g_sigma_coov_ocov_no10_x0_type1_eri_v &
  (s_a, i_a, s_a2, i_a2, T2_, V2_, W34_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_a2, s_a2
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W34_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_a1, i_a1, s_w, i_w, s_a0, i_a0
! W34(w,a0,a1,a) += (    1.00000000) V2(a,a2,v0,a1) T2(w,a0,v0,a2) 
do s_v0 = 0, nir-1
do s_a1 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_a1,s_a) .and. & 
IEOR(s_a,s_a2) == IEOR(s_v0,s_a1) .and. &
IEOR(s_w,s_a0) == IEOR(s_v0,s_a2)) then

if(psym(I_LENGTH,I_O, s_a1) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- V2(a,a2,v0,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z1_(i_a1, i_v0) =  &
  V2_(s_a1, s_v0, s_a2)%array(i_a1, i_v0, i_a2)
end do
end do
! Z2 <-- T2(w,a0,v0,a2) 
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

! Z3 <-- W34(w,a0,a1,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a1),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a1))

! W34(w,a0,a1,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
W34_(s_a1, s_a0, s_w)%array(i_a1, i_a0, i_w) = &
    W34_(s_a1, s_a0, s_w)%array(i_a1, i_a0, i_w) &
  + Z3_(i_a1, i_w, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a1) * &
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

end subroutine g_sigma_coov_ocov_no10_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no10_x1_type1_eri_v &
  (sa, ia, W34, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: W34(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa

call set_symblock_Xcaa(sleft, W34, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no10_x1_type1_eri_v &
  (sa, ia, Xcaa, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no10_x1_type1_eri_v



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
subroutine g_sigma_coov_ocov_no10_x1_type1_eri_v &
  (s_a, i_a, W34_, S2_, D2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W34_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_a1, i_a1, s_a0, i_a0, s_i, i_i, s_w, i_w
! S2(w,i,j,a) += (   -1.00000000) D2(j,a1,a0,i) W34(w,a0,a1,a) 
do s_j = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_j,s_a1) == IEOR(s_a0,s_i) .and. &
IEOR(s_w,s_a0) == IEOR(s_a1,s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(j,a1,a0,i) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a1, i_a0) =  &
  D2_(s_i, s_a0, s_a1, s_j)%array(i_i, i_a0, i_a1, i_j)
end do
end do
end do
end do
! Z2 <-- W34(w,a0,a1,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_w) =  &
  W34_(s_a1, s_a0, s_w)%array(i_a1, i_a0, i_w)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
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

end subroutine g_sigma_coov_ocov_no10_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no11_x0_type1_eri_v &
  (sa, ia, sa1, ia1, T2, V2, W35, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), V2(*), W35(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcaa(sleft, W35, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_coov_ocov_no11_x0_type1_eri_v &
  (sa, ia, sa1, ia1, av2_i, h2_i, Xcaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no11_x0_type1_eri_v



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
subroutine g_sigma_coov_ocov_no11_x0_type1_eri_v &
  (s_a, i_a, s_a1, i_a1, T2_, V2_, W35_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W35_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_a2, i_a2, s_w, i_w, s_a0, i_a0
! W35(w,a0,a2,a) += (    1.00000000) V2(a,v0,a2,a1) T2(w,a0,v0,a1) 
do s_v0 = 0, nir-1
do s_a2 = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_a2,s_a) .and. & 
IEOR(s_a,s_v0) == IEOR(s_a2,s_a1) .and. &
IEOR(s_w,s_a0) == IEOR(s_v0,s_a1)) then

if(psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- V2(a,v0,a2,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z1_(i_a2, i_v0) =  &
  V2_(s_a1, s_a2, s_v0)%array(i_a1, i_a2, i_v0)
end do
end do
! Z2 <-- T2(w,a0,v0,a1) 
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

! Z3 <-- W35(w,a0,a2,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a2),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a2))

! W35(w,a0,a2,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W35_(s_a2, s_a0, s_w)%array(i_a2, i_a0, i_w) = &
    W35_(s_a2, s_a0, s_w)%array(i_a2, i_a0, i_w) &
  + Z3_(i_a2, i_w, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a2) * &
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

end subroutine g_sigma_coov_ocov_no11_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no11_x1_type1_eri_v &
  (sa, ia, W35, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: W35(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa

call set_symblock_Xcaa(sleft, W35, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no11_x1_type1_eri_v &
  (sa, ia, Xcaa, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no11_x1_type1_eri_v



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
subroutine g_sigma_coov_ocov_no11_x1_type1_eri_v &
  (s_a, i_a, W35_, S2_, D2_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W35_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_i, i_i, s_a0, i_a0, s_a2, i_a2, s_w, i_w
! S2(w,i,j,a) += (   -1.00000000) D2(j,i,a0,a2) W35(w,a0,a2,a) 
do s_j = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_a2 = 0, nir-1
do s_w = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_j,s_i) == IEOR(s_a0,s_a2) .and. &
IEOR(s_w,s_a0) == IEOR(s_a2,s_a)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- D2(j,i,a0,a2) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a0, i_a2) =  &
  D2_(s_a2, s_a0, s_i, s_j)%array(i_a2, i_a0, i_i, i_j)
end do
end do
end do
end do
! Z2 <-- W35(w,a0,a2,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
Z2_(i_a0, i_a2, i_w) =  &
  W35_(s_a2, s_a0, s_w)%array(i_a2, i_a0, i_w)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a2),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_j, i_i, i_w)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_C, s_w) * &
                psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_ocov_no11_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no12_x0_type1_eri_v &
  (sa, ia, sa1, ia1, T2, V2, W36, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), V2(*), W36(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcaa(sleft, W36, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_coov_ocov_no12_x0_type1_eri_v &
  (sa, ia, sa1, ia1, av2_i, h2_i, Xcaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no12_x0_type1_eri_v



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
subroutine g_sigma_coov_ocov_no12_x0_type1_eri_v &
  (s_a, i_a, s_a1, i_a1, T2_, V2_, W36_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W36_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_j, i_j, s_w, i_w, s_a0, i_a0
! W36(w,a0,j,a) += (    1.00000000) V2(a,v0,j,a1) T2(w,a0,v0,a1) 
do s_v0 = 0, nir-1
do s_j = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_j,s_a) .and. & 
IEOR(s_a,s_v0) == IEOR(s_j,s_a1) .and. &
IEOR(s_w,s_a0) == IEOR(s_v0,s_a1)) then

if(psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- V2(a,v0,j,a1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_v0) =  &
  V2_(s_a1, s_j, s_v0)%array(i_a1, i_j, i_v0)
end do
end do
! Z2 <-- T2(w,a0,v0,a1) 
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

! Z3 <-- W36(w,a0,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j))

! W36(w,a0,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
W36_(s_j, s_a0, s_w)%array(i_j, i_a0, i_w) = &
    W36_(s_j, s_a0, s_w)%array(i_j, i_a0, i_w) &
  + Z3_(i_j, i_w, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no12_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no12_x1_type1_eri_v &
  (sa, ia, W36, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: W36(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa

call set_symblock_Xcaa(sleft, W36, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no12_x1_type1_eri_v &
  (sa, ia, Xcaa, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no12_x1_type1_eri_v



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
subroutine g_sigma_coov_ocov_no12_x1_type1_eri_v &
  (s_a, i_a, W36_, S2_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W36_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_w, i_w, s_j, i_j
! S2(w,i,j,a) += (   -1.00000000) D1(i,a0) W36(w,a0,j,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_j = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_j,s_a)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W36(w,a0,j,a) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_j, i_a0) =  &
  W36_(s_j, s_a0, s_w)%array(i_j, i_a0, i_w)
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

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     - 1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_w, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no12_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no13_x0_type1_eri_v &
  (sa, ia, sa1, ia1, T2, V2, W37, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sa1, ia1
real(kind=8), intent(inout) :: T2(*), V2(*), W37(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sa1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_Xcaa(sleft, W37, nir, nsym, psym) ! -> Xcaa (allocate) 
call g_sigma_coov_ocov_no13_x0_type1_eri_v &
  (sa, ia, sa1, ia1, av2_i, h2_i, Xcaa, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no13_x0_type1_eri_v



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
subroutine g_sigma_coov_ocov_no13_x0_type1_eri_v &
  (s_a, i_a, s_a1, i_a1, T2_, V2_, W37_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W37_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_v0, i_v0, s_j, i_j, s_w, i_w, s_a0, i_a0
! W37(w,a0,j,a) += (    1.00000000) V2(a,a1,v0,j) T2(w,a0,v0,a1) 
do s_v0 = 0, nir-1
do s_j = 0, nir-1
do s_w = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_w,s_a0) == IEOR(s_j,s_a) .and. & 
IEOR(s_a,s_a1) == IEOR(s_v0,s_j) .and. &
IEOR(s_w,s_a0) == IEOR(s_v0,s_a1)) then

if(psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- V2(a,a1,v0,j) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_v0) =  &
  V2_(s_j, s_v0, s_a1)%array(i_j, i_v0, i_a1)
end do
end do
! Z2 <-- T2(w,a0,v0,a1) 
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

! Z3 <-- W37(w,a0,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_V, s_v0),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j))

! W37(w,a0,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
W37_(s_j, s_a0, s_w)%array(i_j, i_a0, i_w) = &
    W37_(s_j, s_a0, s_w)%array(i_j, i_a0, i_w) &
  + Z3_(i_j, i_w, i_a0)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no13_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_ocov_no13_x1_type1_eri_v &
  (sa, ia, W37, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: W37(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

sleft = sa

call set_symblock_Xcaa(sleft, W37, nir, nsym, psym) ! -> Xcaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_ocov_no13_x1_type1_eri_v &
  (sa, ia, Xcaa, av2_i2, d1, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(Xcaa)

end subroutine g_if_sigma_coov_ocov_no13_x1_type1_eri_v



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
subroutine g_sigma_coov_ocov_no13_x1_type1_eri_v &
  (s_a, i_a, W37_, S2_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W37_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_a0, i_a0, s_w, i_w, s_j, i_j
! S2(w,i,j,a) += (    2.00000000) D1(i,a0) W37(w,a0,j,a) 
do s_i = 0, nir-1
do s_a0 = 0, nir-1
do s_w = 0, nir-1
do s_j = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_i,s_a0) == 0 .and. &
IEOR(s_w,s_a0) == IEOR(s_j,s_a)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- W37(w,a0,j,a) 
allocate(Z1_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
Z1_(i_w, i_j, i_a0) =  &
  W37_(s_j, s_a0, s_w)%array(i_j, i_a0, i_w)
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

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a0),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j))

! S2(w,i,j,a)  <-- Z3
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) = &
    S2_(s_j, s_i, s_w)%array(i_j, i_i, i_w) &
  + Z3_(i_w, i_j, i_i)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) * &
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

end subroutine g_sigma_coov_ocov_no13_x1_type1_eri_v

