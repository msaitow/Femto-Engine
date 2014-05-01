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
subroutine g_if_sigma_g_ccvv_no0_x0_type1_eri_v &
  (sv0, iv0, sv1, iv1, T2, V2, S0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_sigma_g_ccvv_no0_x0_type1_eri_v &
  (sv0, iv0, sv1, iv1, av2_i, h2_i, S0, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_g_ccvv_no0_x0_type1_eri_v



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
subroutine g_sigma_g_ccvv_no0_x0_type1_eri_v &
  (s_v0, i_v0, s_v1, i_v1, T2_, V2_, S0_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_v0, s_v0
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
! S0() += (    4.00000000) T2(c1,c0,v1,v0) V2(v1,c1,v0,c0) 
do s_c1 = 0, nir-1
do s_c0 = 0, nir-1
if( &
IEOR(s_c1,s_c0) == IEOR(s_v1,s_v0) .and. &
IEOR(s_v1,s_c1) == IEOR(s_v0,s_c0)) then

if(psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c1,c0,v1,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z1_(i_c1, i_c0) =  &
  T2_(s_v1, s_c0, s_c1)%array(i_v1, i_c0, i_c1)
end do
end do
! Z2 <-- V2(v1,c1,v0,c0) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z2_(i_c1, i_c0) =  &
  V2_(s_c0, s_v0, s_c1)%array(i_c0, i_v0, i_c1)
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

end subroutine g_sigma_g_ccvv_no0_x0_type1_eri_v



!                      >> makeF90_interface2_new <<                     
! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ccvv_no1_x0_type1_eri_v &
  (sv0, iv0, sv1, iv1, T2, V2, S0, nir, nsym, psym, flops)

use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use f_mr_if_module
implicit none

integer, intent(inout) :: sv0, iv0, sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), S0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft

call set_symblock_av2(sv0, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
call g_sigma_g_ccvv_no1_x0_type1_eri_v &
  (sv0, iv0, sv1, iv1, av2_i, h2_i, S0, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)

end subroutine g_if_sigma_g_ccvv_no1_x0_type1_eri_v



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
subroutine g_sigma_g_ccvv_no1_x0_type1_eri_v &
  (s_v0, i_v0, s_v1, i_v1, T2_, V2_, S0_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_v0, s_v0
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
! S0() += (   -2.00000000) T2(c1,c0,v1,v0) V2(v1,c0,v0,c1) 
do s_c1 = 0, nir-1
do s_c0 = 0, nir-1
if( &
IEOR(s_c1,s_c0) == IEOR(s_v1,s_v0) .and. &
IEOR(s_v1,s_c0) == IEOR(s_v0,s_c1)) then

if(psym(I_LENGTH,I_C, s_c1)*psym(I_LENGTH,I_C, s_c0) > 0) then

! Z1 <-- T2(c1,c0,v1,v0) 
allocate(Z1_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z1_(i_c1, i_c0) =  &
  T2_(s_v1, s_c0, s_c1)%array(i_v1, i_c0, i_c1)
end do
end do
! Z2 <-- V2(v1,c0,v0,c1) 
allocate(Z2_(psym(I_BEGIN,I_C, s_c1):psym(I_END,I_C, s_c1), &
             psym(I_BEGIN,I_C, s_c0):psym(I_END,I_C, s_c0)))
do i_c0 = psym(I_BEGIN, I_C, s_c0), psym(I_END, I_C, s_c0)
do i_c1 = psym(I_BEGIN, I_C, s_c1), psym(I_END, I_C, s_c1)
Z2_(i_c1, i_c0) =  &
  V2_(s_c1, s_v0, s_c0)%array(i_c1, i_v0, i_c0)
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

end subroutine g_sigma_g_ccvv_no1_x0_type1_eri_v

