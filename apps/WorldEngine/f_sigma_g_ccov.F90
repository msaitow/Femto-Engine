#include <sci/icmr/fsrc/f_mr.fh>


!  8888888888                     888                  
!  888                            888                  
!  888                            888                  
!  8888888  .d88b.  88888b.d88b.  888888  .d88b.       
!  888     d8P  Y8b 888 "888 "88b 888    d88""88b  
!  888     88888888 888  888  888 888    888  888      
!  888     Y8b.     888  888  888 Y88b.  Y88..88P      
!  888      "Y8888  888  888  888  "Y888  "Y88P"   

!                                    Generated date : Sun Apr 20 10:26:26 2014



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ccov_no0_x0_type1_eri_v &
  (sv0, iv0, V2, W0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: V2(*), W0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv0

call set_symblock_Xcca(sleft, W0, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_g_ccov_no0_x0_type1_eri_v &
  (sv0, iv0, h2_i, Xcca, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcca)

end subroutine g_if_sigma_g_ccov_no0_x0_type1_eri_v



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
subroutine g_sigma_g_ccov_no0_x0_type1_eri_v &
  (s_v0, i_v0, V2_, W0_, D1_, nir, nsym, psym, flops)

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
type(symblock3), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_c0, i_c0, s_c1, i_c1, s_a1, i_a1, s_a0, i_a0
! W0(c1,c0,a0,v0) += (    1.00000000) V2(v0,c0,c1,a1) D1(a1,a0) 
do s_c0 = 0, nir-1
do s_c1 = 0, nir-1
do s_a1 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_c1,s_c0) == IEOR(s_a0,s_v0) .and. & 
IEOR(s_v0,s_c0) == IEOR(s_c1,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_O, s_a0) > 0 .and. &
   psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1) > 0 .and. &
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
! Z2 <-- V2(v0,c0,c1,a1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1)))
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1, i_c0, i_c1) =  &
  V2_(s_a1, s_c1, s_c0)%array(i_a1, i_c1, i_c0)
end do
end do
end do

! Z3 <-- W0(c1,c0,a0,v0) 
allocate(Z3_(psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_a0),&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1),&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_a0),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_a0))

! W0(c1,c0,a0,v0)  <-- Z3
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
W0_(s_a0, s_c0, s_c1)%array(i_a0, i_c0, i_c1) = &
    W0_(s_a0, s_c0, s_c1)%array(i_a0, i_c0, i_c1) &
  + Z3_(i_a0, i_c0, i_c1)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_a0) * &
                psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1) * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_g_ccov_no0_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ccov_no0_x1_type1_eri_v &
  (sv0, iv0, T2, W0, S0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0
real(kind=8), intent(inout) :: T2(*), W0(*), S0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv0

call set_symblock_Xcca(sleft, W0, nir, nsym, psym) ! -> Xcca (allocate) 
call g_sigma_g_ccov_no0_x1_type1_eri_v &
  (sv0, iv0, av2_i, Xcca, S0, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcca)

end subroutine g_if_sigma_g_ccov_no0_x1_type1_eri_v



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
subroutine g_sigma_g_ccov_no0_x1_type1_eri_v &
  (s_v0, i_v0, T2_, W0_, S0_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
real(kind=8)                   :: S0_
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_c0, i_c0, s_a0, i_a0
! S0() += (   -2.00000000) T2(c1,c0,a0,v0) W0(c1,c0,a0,v0) 
do s_c1 = 0, nir-1
do s_c0 = 0, nir-1
do s_a0 = 0, nir-1
if( &
IEOR(s_c1,s_c0) == IEOR(s_a0,s_v0) .and. &
IEOR(s_c1,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) > 0) then

! Z1 <-- T2(c1,c0,a0,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z1_(i_c1, i_c0, i_a0) =  &
  T2_(s_a0, s_c0, s_c1)%array(i_a0, i_c0, i_c1)
end do
end do
end do
! Z2 <-- W0(c1,c0,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_O, s_a0):psym(I_END,I_O, s_a0)))
do i_a0 = psym(I_BEGIN, I_O, s_a0), psym(I_END, I_O, s_a0)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z2_(i_c1, i_c0, i_a0) =  &
  W0_(s_a0, s_c0, s_c1)%array(i_a0, i_c0, i_c1)
end do
end do
end do

! Gemm Z1 * Z2 to form S0() 
call dgemm('n', 'n', 1,&
                     1,&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     1,&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0),&
                     1.0d+00,&
                     S0_,&
                     1)

! Flop count
flops = flops + 1 * &
                1 * &
                psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_O, s_a0) * 2.0d+00

deallocate(Z1_, Z2_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_g_ccov_no0_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ccov_no1_x0_type1_eri_v &
  (sa0, ia0, sv0, iv0, V2, W1, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0, sv0, iv0
real(kind=8), intent(inout) :: V2(*), W1(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sa0,sv0)

call set_symblock_Xcc(sleft, W1, nir, nsym, psym) ! -> Xcc (allocate) 
call g_sigma_g_ccov_no1_x0_type1_eri_v &
  (sa0, ia0, sv0, iv0, h2_i, Xcc, d1, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xcc)

end subroutine g_if_sigma_g_ccov_no1_x0_type1_eri_v



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
subroutine g_sigma_g_ccov_no1_x0_type1_eri_v &
  (s_a0, i_a0, s_v0, i_v0, V2_, W1_, D1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_a0, s_a0
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

integer :: s_c0, i_c0, s_c1, i_c1, s_a1, i_a1
! W1(c1,c0,a0,v0) += (    1.00000000) V2(v0,c0,c1,a1) D1(a1,a0) 
do s_c0 = 0, nir-1
do s_c1 = 0, nir-1
do s_a1 = 0, nir-1
if( &
IEOR(s_c1,s_c0) == IEOR(s_a0,s_v0) .and. & 
IEOR(s_v0,s_c0) == IEOR(s_c1,s_a1) .and. &
IEOR(s_a1,s_a0) == 0) then

if(psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1) > 0 .and. &
   psym(I_LENGTH,I_O, s_a1) > 0) then

! Z1 <-- V2(v0,c0,c1,a1) 
allocate(Z1_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
Z1_(i_c0, i_c1, i_a1) =  &
  V2_(s_a1, s_c1, s_c0)%array(i_a1, i_c1, i_c0)
end do
end do
end do
! Z2 <-- D1(a1,a0) 
allocate(Z2_(psym(I_BEGIN,I_O, s_a1):psym(I_END,I_O, s_a1)))
do i_a1 = psym(I_BEGIN, I_O, s_a1), psym(I_END, I_O, s_a1)
Z2_(i_a1) =  &
  D1_(s_a0, s_a1)%array(i_a0, i_a1)
end do

! Z3 <-- W1(c1,c0,a0,v0) 
allocate(Z3_(psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0), &
             psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1),&
                     1,&
                     psym(I_LENGTH,I_O, s_a1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_a1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1))

! W1(c1,c0,a0,v0)  <-- Z3
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
W1_(s_c0, s_c1)%array(i_c0, i_c1) = &
    W1_(s_c0, s_c1)%array(i_c0, i_c1) &
  + Z3_(i_c0, i_c1)
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_C, s_c0)*psym(I_LENGTH,I_C, s_c1) * &
                1 * &
                psym(I_LENGTH,I_O, s_a1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_g_ccov_no1_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ccov_no1_x1_type1_eri_v &
  (sa0, ia0, sv0, iv0, T2, W1, S0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0, sv0, iv0
real(kind=8), intent(inout) :: T2(*), W1(*), S0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sa0,sv0)

call set_symblock_Xcc(sleft, W1, nir, nsym, psym) ! -> Xcc (allocate) 
call g_sigma_g_ccov_no1_x1_type1_eri_v &
  (sa0, ia0, sv0, iv0, av2_i, Xcc, S0, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xcc)

end subroutine g_if_sigma_g_ccov_no1_x1_type1_eri_v



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
subroutine g_sigma_g_ccov_no1_x1_type1_eri_v &
  (s_a0, i_a0, s_v0, i_v0, T2_, W1_, S0_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
real(kind=8)                   :: S0_
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: W1_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_c0, i_c0
! S0() += (    1.00000000) T2(c1,c0,v0,a0) W1(c1,c0,a0,v0) 
do s_c1 = 0, nir-1
do s_c0 = 0, nir-1
if( &
IEOR(s_c1,s_c0) == IEOR(s_v0,s_a0) .and. &
IEOR(s_c1,s_c0) == IEOR(s_a0,s_v0)) then

if(psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c1,c0,v0,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z1_(i_c1, i_c0) =  &
  T2_(s_v0, s_c0, s_c1)%array(i_v0, i_c0, i_c1)
end do
end do
! Z2 <-- W1(c1,c0,a0,v0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z2_(i_c1, i_c0) =  &
  W1_(s_c0, s_c1)%array(i_c0, i_c1)
end do
end do

! Gemm Z1 * Z2 to form S0() 
call dgemm('n', 'n', 1,&
                     1,&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0),&
                     1.00000000d+00, &
                     Z1_,&
                     1,&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0),&
                     1.0d+00,&
                     S0_,&
                     1)

! Flop count
flops = flops + 1 * &
                1 * &
                psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_g_ccov_no1_x1_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ccov_no2_x0_type1_eri_v &
  (sa0, ia0, sv0, iv0, T2, V2, S0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0, sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), S0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_sigma_g_ccov_no2_x0_type1_eri_v &
  (sa0, ia0, sv0, iv0, av2_i, h2_i, S0, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_g_ccov_no2_x0_type1_eri_v



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
subroutine g_sigma_g_ccov_no2_x0_type1_eri_v &
  (s_a0, i_a0, s_v0, i_v0, T2_, V2_, S0_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
real(kind=8)                   :: S0_
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_c0, i_c0
! S0() += (    4.00000000) T2(c1,c0,v0,a0) V2(v0,c1,c0,a0) 
do s_c1 = 0, nir-1
do s_c0 = 0, nir-1
if( &
IEOR(s_c1,s_c0) == IEOR(s_v0,s_a0) .and. &
IEOR(s_v0,s_c1) == IEOR(s_c0,s_a0)) then

if(psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c1,c0,v0,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z1_(i_c1, i_c0) =  &
  T2_(s_v0, s_c0, s_c1)%array(i_v0, i_c0, i_c1)
end do
end do
! Z2 <-- V2(v0,c1,c0,a0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z2_(i_c1, i_c0) =  &
  V2_(s_a0, s_c0, s_c1)%array(i_a0, i_c0, i_c1)
end do
end do

! Gemm Z1 * Z2 to form S0() 
call dgemm('n', 'n', 1,&
                     1,&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0),&
                     4.00000000d+00, &
                     Z1_,&
                     1,&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0),&
                     1.0d+00,&
                     S0_,&
                     1)

! Flop count
flops = flops + 1 * &
                1 * &
                psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_g_ccov_no2_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ccov_no3_x0_type1_eri_v &
  (sa0, ia0, sv0, iv0, T2, V2, S0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sa0, ia0, sv0, iv0
real(kind=8), intent(inout) :: T2(*), V2(*), S0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sa0, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv0, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_sigma_g_ccov_no3_x0_type1_eri_v &
  (sa0, ia0, sv0, iv0, av2_i, h2_i, S0, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_g_ccov_no3_x0_type1_eri_v



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
subroutine g_sigma_g_ccov_no3_x0_type1_eri_v &
  (s_a0, i_a0, s_v0, i_v0, T2_, V2_, S0_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v0, s_v0
integer, intent(in) :: i_a0, s_a0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
real(kind=8)                   :: S0_
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:)
real*8 :: Z3_
! Indices used in the contractions as dummy ... 

integer :: s_c1, i_c1, s_c0, i_c0
! S0() += (   -2.00000000) T2(c1,c0,v0,a0) V2(v0,c0,c1,a0) 
do s_c1 = 0, nir-1
do s_c0 = 0, nir-1
if( &
IEOR(s_c1,s_c0) == IEOR(s_v0,s_a0) .and. &
IEOR(s_v0,s_c0) == IEOR(s_c1,s_a0)) then

if(psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c1,c0,v0,a0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z1_(i_c1, i_c0) =  &
  T2_(s_v0, s_c0, s_c1)%array(i_v0, i_c0, i_c1)
end do
end do
! Z2 <-- V2(v0,c0,c1,a0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z2_(i_c1, i_c0) =  &
  V2_(s_a0, s_c1, s_c0)%array(i_a0, i_c1, i_c0)
end do
end do

! Gemm Z1 * Z2 to form S0() 
call dgemm('n', 'n', 1,&
                     1,&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0),&
                     - 2.00000000d+00, &
                     Z1_,&
                     1,&
                     Z2_,&
                     psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0),&
                     1.0d+00,&
                     S0_,&
                     1)

! Flop count
flops = flops + 1 * &
                1 * &
                psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0) * 2.0d+00

deallocate(Z1_, Z2_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_g_ccov_no3_x0_type1_eri_v

