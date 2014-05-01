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
subroutine g_if_sigma_coov_oovv_no0_x0_type1_eri_c &
  (sw, iw, V2, W0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sw, iw
real(kind=8), intent(inout) :: V2(*), W0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sw, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sw

call set_symblock_Xaaaav(sleft, W0, nir, nsym, psym) ! -> Xaaaav (allocate) 
call g_sigma_coov_oovv_no0_x0_type1_eri_c &
  (sw, iw, h2_i, Xaaaav, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaaaav)

end subroutine g_if_sigma_coov_oovv_no0_x0_type1_eri_c



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
! >> Score :: 110110
! RDM is rotated :: D3(j,a3,a1,i,a0,a2)  >> D3(i,a1,a3,j,a2,a0) 
! summedInd : @[a3, "active"] @[a2, "active"] 
! colInd : @[j, "active"] @[a0, "active"] @[a1, "active"] @[i, "active"] 
subroutine g_sigma_coov_oovv_no0_x0_type1_eri_c &
  (s_w, i_w, V2_, W0_, D3_, nir, nsym, psym, flops)

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
type(symblock5), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:,:,:)
real*8, allocatable :: Z3_(:,:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a3, i_a3, s_v0, i_v0, s_a2, i_a2, s_j, i_j, s_a1, i_a1
integer :: s_i, i_i, s_a0, i_a0
! W0(w,j,a1,a0,i,v0) += (    1.00000000) V2(w,a3,v0,a2) D3(j,a3,a1,i,a0,a2) 
do s_a3 = 0, nir-1
do s_v0 = 0, nir-1
do s_a2 = 0, nir-1
do s_j = 0, nir-1
do s_a1 = 0, nir-1
do s_i = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(IEOR(s_w,s_j),s_a1) == IEOR(IEOR(s_a0,s_i),s_v0) .and. & 
IEOR(s_w,s_a3) == IEOR(s_v0,s_a2) .and. &
IEOR(IEOR(s_j,s_a3),s_a1) == IEOR(IEOR(s_i,s_a0),s_a2)) then

if(psym(I_LENGTH,I_V, s_v0) > 0 .and. &
   psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(w,a3,v0,a2) 
allocate(Z1_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
Z1_(i_v0, i_a3, i_a2) =  &
  V2_(s_a2, s_v0, s_a3)%array(i_a2, i_v0, i_a3)
end do
end do
end do
! Z2 <-- D3(i,a1,a3,j,a2,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a3):psym(I_END,I_O, s_a3), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_a3 = psym(I_BEGIN, I_O, s_a3), psym(I_END, I_O, s_a3)
Z2_(i_a3, i_a2, i_j, i_a0, i_a1, i_i) =  &
  D3_(s_a0, s_a2, s_j, s_a3, s_a1, s_i)%array(i_a0, i_a2, i_j, i_a3, i_a1, i_i)
end do
end do
end do
end do
end do
end do

! Z3 <-- W0(w,j,a1,a0,i,v0) 
allocate(Z3_(psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0), &
             psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_v0),&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_v0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_v0))

! W0(w,j,a1,a0,i,v0)  <-- Z3
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
W0_(s_v0, s_i, s_a0, s_a1, s_j)%array(i_v0, i_i, i_a0, i_a1, i_j) = &
    W0_(s_v0, s_i, s_a0, s_a1, s_j)%array(i_v0, i_i, i_a0, i_a1, i_j) &
  + Z3_(i_v0, i_j, i_a0, i_a1, i_i)
end do
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_v0) * &
                psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_i) * &
                psym(I_LENGTH,I_O, s_a3)*psym(I_LENGTH,I_O, s_a2) * 2.0d+00

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

end subroutine g_sigma_coov_oovv_no0_x0_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_oovv_no0_x1_type1_eri_c &
  (sa, ia, sw, iw, T2, W0, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sw, iw
real(kind=8), intent(inout) :: T2(*), W0(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sw

call set_symblock_Xaaaav(sleft, W0, nir, nsym, psym) ! -> Xaaaav (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_oovv_no0_x1_type1_eri_c &
  (sa, ia, sw, iw, av2_i, Xaaaav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(av2_i)
deallocate(Xaaaav)

end subroutine g_if_sigma_coov_oovv_no0_x1_type1_eri_c



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
subroutine g_sigma_coov_oovv_no0_x1_type1_eri_c &
  (s_a, i_a, s_w, i_w, T2_, W0_, S2_, nir, nsym, psym, flops)

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
type(symblock5), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a0, i_a0, s_a1, i_a1, s_v0, i_v0, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (   -2.00000000) T2(a0,a1,v0,a) W0(w,j,a1,a0,i,v0) 
do s_a0 = 0, nir-1
do s_a1 = 0, nir-1
do s_v0 = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_a0,s_a1) == IEOR(s_v0,s_a) .and. &
IEOR(IEOR(s_w,s_j),s_a1) == IEOR(IEOR(s_a0,s_i),s_v0)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) > 0) then

! Z1 <-- W0(w,j,a1,a0,i,v0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_i, i_a1, i_a0, i_v0) =  &
  W0_(s_v0, s_i, s_a0, s_a1, s_j)%array(i_v0, i_i, i_a0, i_a1, i_j)
end do
end do
end do
end do
end do
! Z2 <-- T2(a0,a1,v0,a) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_V, s_v0):psym(I_END,I_V, s_v0)))
do i_v0 = psym(I_BEGIN, I_V, s_v0), psym(I_END, I_V, s_v0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0, i_v0) =  &
  T2_(s_v0, s_a1, s_a0)%array(i_v0, i_a1, i_a0)
end do
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_O, s_i),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0),&
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
                psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0)*psym(I_LENGTH,I_V, s_v0) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_coov_oovv_no0_x1_type1_eri_c



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_oovv_no0_x0_type1_eri_v &
  (sa, ia, sv0, iv0, T2, W1, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W1(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa,sv0)

call set_symblock_Xaa(sleft, W1, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_coov_oovv_no0_x0_type1_eri_v &
  (sa, ia, sv0, iv0, av2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_oovv_no0_x0_type1_eri_v



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
subroutine g_sigma_coov_oovv_no0_x0_type1_eri_v &
  (s_a, i_a, s_v0, i_v0, T2_, W1_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W1_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a0, i_a0, s_i, i_i, s_a2, i_a2
! W1(i,a2,a,v0) += (    1.00000000) T2(a1,a0,a,v0) D2(i,a1,a2,a0) 
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_a2 = 0, nir-1
if( &
IEOR(s_i,s_a2) == IEOR(s_a,s_v0) .and. & 
IEOR(s_a1,s_a0) == IEOR(s_a,s_v0) .and. &
IEOR(s_i,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(i,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a2, i_a1, i_a0) =  &
  D2_(s_a0, s_a2, s_a1, s_i)%array(i_a0, i_a2, i_a1, i_i)
end do
end do
end do
end do
! Z2 <-- T2(a1,a0,a,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  T2_(s_a, s_a0, s_a1)%array(i_a, i_a0, i_a1)
end do
end do

! Z3 <-- W1(i,a2,a,v0) 
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

! W1(i,a2,a,v0)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W1_(s_a2, s_i)%array(i_a2, i_i) = &
    W1_(s_a2, s_i)%array(i_a2, i_i) &
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

end subroutine g_sigma_coov_oovv_no0_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_oovv_no0_x1_type1_eri_v &
  (sa, ia, sv0, iv0, V2, W1, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W1(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa,sv0)

call set_symblock_Xaa(sleft, W1, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_oovv_no0_x1_type1_eri_v &
  (sa, ia, sv0, iv0, h2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_oovv_no0_x1_type1_eri_v



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
subroutine g_sigma_coov_oovv_no0_x1_type1_eri_v &
  (s_a, i_a, s_v0, i_v0, V2_, W1_, S2_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W1_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a2, i_a2, s_w, i_w, s_j, i_j, s_i, i_i
! S2(w,i,j,a) += (    4.00000000) V2(v0,a2,w,j) W1(i,a2,a,v0) 
do s_a2 = 0, nir-1
do s_w = 0, nir-1
do s_j = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_v0,s_a2) == IEOR(s_w,s_j) .and. &
IEOR(s_i,s_a2) == IEOR(s_a,s_v0)) then

if(psym(I_LENGTH,I_C, s_w)*psym(I_LENGTH,I_O, s_j) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(v0,a2,w,j) 
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
! Z2 <-- W1(i,a2,a,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_i) =  &
  W1_(s_a2, s_i)%array(i_a2, i_i)
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
                     4.00000000d+00, &
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

end subroutine g_sigma_coov_oovv_no0_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_oovv_no1_x0_type1_eri_v &
  (sa, ia, sv0, iv0, T2, W2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa,sv0)

call set_symblock_Xaa(sleft, W2, nir, nsym, psym) ! -> Xaa (allocate) 
call g_sigma_coov_oovv_no1_x0_type1_eri_v &
  (sa, ia, sv0, iv0, av2_i, Xaa, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_oovv_no1_x0_type1_eri_v



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
subroutine g_sigma_coov_oovv_no1_x0_type1_eri_v &
  (s_a, i_a, s_v0, i_v0, T2_, W2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W2_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a1, i_a1, s_a0, i_a0, s_i, i_i, s_a2, i_a2
! W2(i,a2,a,v0) += (    1.00000000) T2(a1,a0,a,v0) D2(i,a1,a2,a0) 
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
do s_i = 0, nir-1
do s_a2 = 0, nir-1
if( &
IEOR(s_i,s_a2) == IEOR(s_a,s_v0) .and. & 
IEOR(s_a1,s_a0) == IEOR(s_a,s_v0) .and. &
IEOR(s_i,s_a1) == IEOR(s_a2,s_a0)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_a2) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- D2(i,a1,a2,a0) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
Z1_(i_i, i_a2, i_a1, i_a0) =  &
  D2_(s_a0, s_a2, s_a1, s_i)%array(i_a0, i_a2, i_a1, i_i)
end do
end do
end do
end do
! Z2 <-- T2(a1,a0,a,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_a0) =  &
  T2_(s_a, s_a0, s_a1)%array(i_a, i_a0, i_a1)
end do
end do

! Z3 <-- W2(i,a2,a,v0) 
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

! W2(i,a2,a,v0)  <-- Z3
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
W2_(s_a2, s_i)%array(i_a2, i_i) = &
    W2_(s_a2, s_i)%array(i_a2, i_i) &
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

end subroutine g_sigma_coov_oovv_no1_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_coov_oovv_no1_x1_type1_eri_v &
  (sa, ia, sv0, iv0, V2, W2, S2, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa, ia, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa,sv0)

call set_symblock_Xaa(sleft, W2, nir, nsym, psym) ! -> Xaa (allocate) 
call set_symblock_av2_2(sa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_coov_oovv_no1_x1_type1_eri_v &
  (sa, ia, sv0, iv0, h2_i, Xaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i2)
deallocate(h2_i)
deallocate(Xaa)

end subroutine g_if_sigma_coov_oovv_no1_x1_type1_eri_v



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
subroutine g_sigma_coov_oovv_no1_x1_type1_eri_v &
  (s_a, i_a, s_v0, i_v0, V2_, W2_, S2_, nir, nsym, psym, flops)

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
type(symblock2), intent(inout) :: W2_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_j, i_j, s_w, i_w, s_a2, i_a2, s_i, i_i
! S2(w,i,j,a) += (   -2.00000000) V2(v0,j,w,a2) W2(i,a2,a,v0) 
do s_j = 0, nir-1
do s_w = 0, nir-1
do s_a2 = 0, nir-1
do s_i = 0, nir-1
if( &
IEOR(s_w,s_i) == IEOR(s_j,s_a) .and. & 
IEOR(s_v0,s_j) == IEOR(s_w,s_a2) .and. &
IEOR(s_i,s_a2) == IEOR(s_a,s_v0)) then

if(psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_C, s_w) > 0 .and. &
   psym(I_LENGTH,I_O, s_i) > 0 .and. &
   psym(I_LENGTH,I_O, s_a2) > 0) then

! Z1 <-- V2(v0,j,w,a2) 
allocate(Z1_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2)))
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
do i_w = psym(I_BEGIN, I_C, s_w), psym(I_END, I_C, s_w)
do i_j = psym(I_BEGIN, I_O, s_j), psym(I_END, I_O, s_j)
Z1_(i_j, i_w, i_a2) =  &
  V2_(s_a2, s_w, s_j)%array(i_a2, i_w, i_j)
end do
end do
end do
! Z2 <-- W2(i,a2,a,v0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a2):psym(I_END,I_O, s_a2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_a2 = psym(I_BEGIN, I_O, s_a2), psym(I_END, I_O, s_a2)
Z2_(i_a2, i_i) =  &
  W2_(s_a2, s_i)%array(i_a2, i_i)
end do
end do

! Z3 <-- S2(w,i,j,a) 
allocate(Z3_(psym(I_BEGIN,I_O, s_j):psym(I_END,I_O, s_j), &
             psym(I_BEGIN,I_C, s_w):psym(I_END,I_C, s_w), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_C, s_w),&
                     psym(I_LENGTH,I_O, s_i),&
                     psym(I_LENGTH,I_O, s_a2),&
                     - 2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_C, s_w),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_C, s_w))

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
flops = flops + psym(I_LENGTH,I_O, s_j)*psym(I_LENGTH,I_C, s_w) * &
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

end subroutine g_sigma_coov_oovv_no1_x1_type1_eri_v

