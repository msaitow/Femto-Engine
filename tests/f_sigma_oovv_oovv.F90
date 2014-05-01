#include "../f_ct.fh"


!  ___________                __               
!  \_   _____/____    _____ _/  |_  ____      
!   |    __)_/ __ \  /     \\   __\/  _ \ 
!   |     \ \  ___/ |  Y Y  \|  | (  <_> )  
!   \___  /  \___  >|__|_|  /|__|  \____/   
!       \/       \/       \/                



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x0_type1_noeri(sc, ic, Fc0, T2, S2, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: Fc0
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no0_x0_type1_noeri(sc, ic, Fc0, av2_i, av2_i2, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no0_x0_type1_noeri



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x0_type1_noeri(s_c, i_c, Fc0, T2_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Fc0
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2, s_a, i_a
! S2(i,k,a,c) += (    2.00000000) Fc0 D2(i,o1,k,o2) T2(o1,o2,a,c) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_o2) .and. &
IEOR(s_o1,s_o2) == IEOR(s_a,s_c)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2) > 0) then

! Z1 <-- T2(o1,o2,a,c) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2)))
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Z1_(i_a, i_o1, i_o2) =  &
  T2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
end do
end do
end do
! Z2 <-- D2(i,o1,k,o2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_o1, i_o2, i_i, i_k) =  &
  D2_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i)
end do
end do
end do
end do

! Z3 <-- S2(i,k,a,c) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2),&
                     2.00000000d+00*Fc0, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,k,a,c)  <-- Z3
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) = &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + Z3_(i_a, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no0_x0_type1_noeri



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x1_type1_noeri(W0, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: W0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = 0
call set_symblock_Xaaaa(sleft, W0, nir, nsym, psym) ! -> Xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x1_type1_noeri(Xaaaa, d3, fc1, nir, nsym, psym, flops)

deallocate(Xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x1_type1_noeri



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x1_type1_noeri(W0_, D3_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2, s_o3, i_o3
integer :: s_o4, i_o4
! W0(i,o1,k,o2) += (    1.00000000) D3(i,o1,k,o2,o3,o4) Fc1(o3,o4) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_i,s_o1) == IEOR(s_k,s_o2) .and. & 
IEOR(IEOR(s_i,s_o1),s_k) == IEOR(IEOR(s_o2,s_o3),s_o4) .and. &
IEOR(s_o3,s_o4) == 0) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o2) > 0 .and. &
   psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4) > 0) then

! Z1 <-- D3(i,o1,k,o2,o3,o4) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Z1_(i_i, i_o1, i_k, i_o2, i_o3, i_o4) =  &
  D3_(s_o4, s_o3, s_o2, s_k, s_o1, s_i)%array(i_o4, i_o3, i_o2, i_k, i_o1, i_i)
end do
end do
end do
end do
end do
end do
! Z2 <-- Fc1(o3,o4) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4)))
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Z2_(i_o3, i_o4) =  &
  Fc1_(s_o4, s_o3)%array(i_o4, i_o3)
end do
end do

! Z3 <-- W0(i,o1,k,o2) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o2),&
                     1,&
                     psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o2))

! W0(i,o1,k,o2)  <-- Z3
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W0_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i) = &
    W0_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i) &
  + Z3_(i_i, i_o1, i_k, i_o2)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o2) * &
                1 * &
                psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4) * 2.0d+00

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

end subroutine g_sigma_oovv_oovv_no0_x1_type1_noeri



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x1_type1_noeri(sc, ic, T2, W0, S2, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), W0(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_Xaaaa(sleft, W0, nir, nsym, psym) ! -> Xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x1_type1_noeri(sc, ic, av2_i, Xaaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x1_type1_noeri



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x1_type1_noeri(s_c, i_c, T2_, W0_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_a, i_a, s_i, i_i, s_k, i_k
! S2(i,k,a,c) += (    2.00000000) T2(o1,o2,a,c) W0(i,o1,k,o2) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_a,s_c) .and. &
IEOR(s_i,s_o1) == IEOR(s_k,s_o2)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2) > 0) then

! Z1 <-- T2(o1,o2,a,c) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2)))
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Z1_(i_a, i_o1, i_o2) =  &
  T2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
end do
end do
end do
! Z2 <-- W0(i,o1,k,o2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_o1, i_o2, i_i, i_k) =  &
  W0_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i)
end do
end do
end do
end do

! Z3 <-- S2(i,k,a,c) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,k,a,c)  <-- Z3
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) = &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + Z3_(i_a, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no1_x1_type1_noeri



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x2_type1_noeri(sv1, iv1, T2, W1, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv1, iv1
real(kind=8), intent(inout) :: T2(*), W1(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv1

call set_symblock_Xaav(sleft, W1, nir, nsym, psym) ! -> Xaav (allocate) 
call g_sigma_oovv_oovv_no0_x2_type1_noeri(sv1, iv1, av2_i, Xaav, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x2_type1_noeri



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x2_type1_noeri(s_v1, i_v1, T2_, W1_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: W1_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2, s_a, i_a
! W1(i,k,a,v1) += (    1.00000000) D2(i,o1,k,o2) T2(o1,o2,a,v1) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_v1) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_o2) .and. &
IEOR(s_o1,s_o2) == IEOR(s_a,s_v1)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2) > 0) then

! Z1 <-- T2(o1,o2,a,v1) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2)))
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Z1_(i_a, i_o1, i_o2) =  &
  T2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
end do
end do
end do
! Z2 <-- D2(i,o1,k,o2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_o1, i_o2, i_i, i_k) =  &
  D2_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i)
end do
end do
end do
end do

! Z3 <-- W1(i,k,a,v1) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W1(i,k,a,v1)  <-- Z3
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W1_(s_a, s_k, s_i)%array(i_a, i_k, i_i) = &
    W1_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + Z3_(i_a, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no0_x2_type1_noeri



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x2_type1_noeri(sc, ic, sv1, iv1, W1, S2, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: W1(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sv1

call set_symblock_Xaav(sleft, W1, nir, nsym, psym) ! -> Xaav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x2_type1_noeri(sc, ic, sv1, iv1, Xaav, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(Xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x2_type1_noeri



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x2_type1_noeri(s_c, i_c, s_v1, i_v1, W1_, S2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W1_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8 :: Z2_
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_a, i_a
! S2(i,k,a,c) += (    2.00000000) Fc1(c,v1) W1(i,k,a,v1) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_c,s_v1) == 0 .and. &
IEOR(s_i,s_k) == IEOR(s_a,s_v1)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a) > 0) then

! Z1 <-- W1(i,k,a,v1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Z1_(i_i, i_k, i_a) =  &
  W1_(s_a, s_k, s_i)%array(i_a, i_k, i_i)
end do
end do
end do
! Z2 <-- Fc1(c,v1) 
Z2_ =  &
  Fc1_(s_v1, s_c)%array(i_v1, i_c)

! Z3 <-- S2(i,k,a,c) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a),&
                     1,&
                     1,&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     1,&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a))

! S2(i,k,a,c)  <-- Z3
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) = &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + Z3_(i_i, i_k, i_a)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a) * &
                1 * &
                1 * 2.0d+00

deallocate(Z1_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no1_x2_type1_noeri



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x3_type1_noeri(sc, ic, T2, W2, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), W2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc

call set_symblock_Xaav(sleft, W2, nir, nsym, psym) ! -> Xaav (allocate) 
call g_sigma_oovv_oovv_no0_x3_type1_noeri(sc, ic, av2_i, Xaav, d2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x3_type1_noeri



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x3_type1_noeri(s_c, i_c, T2_, W2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: W2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2, s_v1, i_v1
! W2(i,k,v1,c) += (    1.00000000) D2(i,o1,k,o2) T2(o1,o2,v1,c) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_v1,s_c) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_o2) .and. &
IEOR(s_o1,s_o2) == IEOR(s_v1,s_c)) then

if(psym(I_LENGTH,I_V, s_v1) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2) > 0) then

! Z1 <-- T2(o1,o2,v1,c) 
allocate(Z1_(psym(I_BEGIN,I_V, s_v1):psym(I_END,I_V, s_v1), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2)))
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Z1_(i_v1, i_o1, i_o2) =  &
  T2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1)
end do
end do
end do
! Z2 <-- D2(i,o1,k,o2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_o1, i_o2, i_i, i_k) =  &
  D2_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i)
end do
end do
end do
end do

! Z3 <-- W2(i,k,v1,c) 
allocate(Z3_(psym(I_BEGIN,I_V, s_v1):psym(I_END,I_V, s_v1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_v1),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_v1),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_v1))

! W2(i,k,v1,c)  <-- Z3
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W2_(s_v1, s_k, s_i)%array(i_v1, i_k, i_i) = &
    W2_(s_v1, s_k, s_i)%array(i_v1, i_k, i_i) &
  + Z3_(i_v1, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_v1) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no0_x3_type1_noeri



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x3_type1_noeri(sc, ic, W2, S2, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: W2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sc

call set_symblock_Xaav(sleft, W2, nir, nsym, psym) ! -> Xaav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x3_type1_noeri(sc, ic, Xaav, av2_i2, fc1, nir, nsym, psym, flops)

deallocate(Xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x3_type1_noeri



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x3_type1_noeri(s_c, i_c, W2_, S2_, Fc1_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Fc1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_v1, i_v1, s_i, i_i, s_k, i_k
! S2(i,k,a,c) += (    2.00000000) Fc1(a,v1) W2(i,k,v1,c) 
do s_a = 0, nir-1
do s_v1 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_a,s_v1) == 0 .and. &
IEOR(s_i,s_k) == IEOR(s_v1,s_c)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_V, s_v1) > 0) then

! Z1 <-- Fc1(a,v1) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v1):psym(I_END,I_V, s_v1)))
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Z1_(i_a, i_v1) =  &
  Fc1_(s_v1, s_a)%array(i_v1, i_a)
end do
end do
! Z2 <-- W2(i,k,v1,c) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v1):psym(I_END,I_V, s_v1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_v1, i_i, i_k) =  &
  W2_(s_v1, s_k, s_i)%array(i_v1, i_k, i_i)
end do
end do
end do

! Z3 <-- S2(i,k,a,c) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_V, s_v1),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,k,a,c)  <-- Z3
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) = &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + Z3_(i_a, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_V, s_v1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no1_x3_type1_noeri



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x0_type1_eri_o(so1, io1, so5, io5, V2, W0, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: so1, io1, so5, io5
real(kind=8), intent(inout) :: V2(*), W0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(so1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, W0, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x0_type1_eri_o(so1, io1, so5, io5, h2_i, xaaaa, d4_ij, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x0_type1_eri_o



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x0_type1_eri_o(s_o1, i_o1, s_o5, i_o5, V2_, W0_, D4_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_o1, s_o1
integer, intent(in) :: i_o5, s_o5
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o4, i_o4, s_o2, i_o2
integer :: s_o6, i_o6
! W0(i,o3,k,o4) += (    1.00000000) D4(o1,o5,i,o3,k,o4,o2,o6) V2(o1,o5,o2,o6) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
do s_o6 = 0, nir-1
if( &
IEOR(s_i,s_o3) == IEOR(s_k,s_o4) .and. & 
IEOR(IEOR(s_o1,s_o5),IEOR(s_i,s_o3)) == IEOR(IEOR(s_k,s_o4),IEOR(s_o2,s_o6)) .and. &
IEOR(s_o1,s_o5) == IEOR(s_o2,s_o6)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o4) > 0 .and. &
   psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_o6) > 0) then

! Z1 <-- D4(o1,o5,i,o3,k,o4,o2,o6) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_o6):psym(I_END,I_O, s_o6)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o6 = psym(I_BEGIN, I_O, s_o6), psym(I_END, I_O, s_o6)
Z1_(i_i, i_o3, i_k, i_o4, i_o2, i_o6) =  &
  D4_(s_o6, s_o2, s_o4, s_k, s_o3, s_i)%array(i_o6, i_o2, i_o4, i_k, i_o3, i_i)
end do
end do
end do
end do
end do
end do
! Z2 <-- V2(o1,o5,o2,o6) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_o6):psym(I_END,I_O, s_o6)))
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o6 = psym(I_BEGIN, I_O, s_o6), psym(I_END, I_O, s_o6)
Z2_(i_o2, i_o6) =  &
  V2_(s_o6, s_o2, s_o5)%array(i_o6, i_o2, i_o5)
end do
end do

! Z3 <-- W0(i,o3,k,o4) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o4),&
                     1,&
                     psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_o6),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o4),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_o6),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o4))

! W0(i,o3,k,o4)  <-- Z3
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W0_(s_o4, s_k, s_o3, s_i)%array(i_o4, i_k, i_o3, i_i) = &
    W0_(s_o4, s_k, s_o3, s_i)%array(i_o4, i_k, i_o3, i_i) &
  + Z3_(i_i, i_o3, i_k, i_o4)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o4) * &
                1 * &
                psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_o6) * 2.0d+00

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

end subroutine g_sigma_oovv_oovv_no0_x0_type1_eri_o



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x0_type2_eri_o(sc, ic, T2, W0, S2, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T2(*), W0(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, W0, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no0_x0_type2_eri_o(sc, ic, av2_i, xaaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no0_x0_type2_eri_o



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x0_type2_eri_o(s_c, i_c, T2_, W0_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o4, i_o4, s_a, i_a, s_i, i_i, s_k, i_k
! S2(i,k,a,c) += (    1.00000000) T2(o3,o4,a,c) W0(i,o3,k,o4) 
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o3,s_o4) == IEOR(s_a,s_c) .and. &
IEOR(s_i,s_o3) == IEOR(s_k,s_o4)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4) > 0) then

! Z1 <-- T2(o3,o4,a,c) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4)))
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Z1_(i_a, i_o3, i_o4) =  &
  T2_(s_a, s_o4, s_o3)%array(i_a, i_o4, i_o3)
end do
end do
end do
! Z2 <-- W0(i,o3,k,o4) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_o3, i_o4, i_i, i_k) =  &
  W0_(s_o4, s_k, s_o3, s_i)%array(i_o4, i_k, i_o3, i_i)
end do
end do
end do
end do

! Z3 <-- S2(i,k,a,c) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,k,a,c)  <-- Z3
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) = &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + Z3_(i_a, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no0_x0_type2_eri_o



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x0_type1_eri_v(sc, ic, sv1, iv1, V2, W0, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: V2(*), W0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sc,sv1)

call set_symblock_Xaaaa(sleft, W0, nir, nsym, psym) ! -> Xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x0_type1_eri_v(sc, ic, sv1, iv1, h2_i, Xaaaa, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x0_type1_eri_v



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x0_type1_eri_v(s_c, i_c, s_v1, i_v1, V2_, W0_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! W0(i,o3,k,o2,c,v1) += (    1.00000000) D3(i,o3,k,o2,o4,o1) V2(c,v1,o1,o4) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_c),s_v1) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(s_c,s_v1) == IEOR(s_o1,s_o4)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o2) > 0 .and. &
   psym(I_LENGTH,I_O, s_o4)*psym(I_LENGTH,I_O, s_o1) > 0) then

! Z1 <-- D3(i,o3,k,o2,o4,o1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Z1_(i_i, i_o3, i_k, i_o2, i_o4, i_o1) =  &
  D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i)
end do
end do
end do
end do
end do
end do
! Z2 <-- V2(c,v1,o1,o4) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1)))
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Z2_(i_o4, i_o1) =  &
  V2_(s_o4, s_o1, s_v1)%array(i_o4, i_o1, i_v1)
end do
end do

! Z3 <-- W0(i,o3,k,o2,c,v1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o2),&
                     1,&
                     psym(I_LENGTH,I_O, s_o4)*psym(I_LENGTH,I_O, s_o1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o2),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o4)*psym(I_LENGTH,I_O, s_o1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o2))

! W0(i,o3,k,o2,c,v1)  <-- Z3
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W0_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) = &
    W0_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i) &
  + Z3_(i_i, i_o3, i_k, i_o2)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o2) * &
                1 * &
                psym(I_LENGTH,I_O, s_o4)*psym(I_LENGTH,I_O, s_o1) * 2.0d+00

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

end subroutine g_sigma_oovv_oovv_no0_x0_type1_eri_v



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x0_type1_eri_v(sc, ic, sv1, iv1, T2, W0, S2, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), W0(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sc,sv1)

call set_symblock_Xaaaa(sleft, W0, nir, nsym, psym) ! -> Xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x0_type1_eri_v(sc, ic, sv1, iv1, av2_i, Xaaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x0_type1_eri_v



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x0_type1_eri_v(s_c, i_c, s_v1, i_v1, T2_, W0_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W0_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_a, i_a, s_i, i_i, s_k, i_k
! S2(i,k,a,c) += (    2.00000000) T2(o3,o2,a,v1) W0(i,o3,k,o2,c,v1) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_a,s_v1) .and. &
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_c),s_v1)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o2) > 0) then

! Z1 <-- T2(o3,o2,a,v1) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2)))
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Z1_(i_a, i_o3, i_o2) =  &
  T2_(s_a, s_o2, s_o3)%array(i_a, i_o2, i_o3)
end do
end do
end do
! Z2 <-- W0(i,o3,k,o2,c,v1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_o3, i_o2, i_i, i_k) =  &
  W0_(s_o2, s_k, s_o3, s_i)%array(i_o2, i_k, i_o3, i_i)
end do
end do
end do
end do

! Z3 <-- S2(i,k,a,c) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o2),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,k,a,c)  <-- Z3
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) = &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + Z3_(i_a, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no1_x0_type1_eri_v



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x1_type1_eri_v(sv1, iv1, V2, W1, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv1, iv1
real(kind=8), intent(inout) :: V2(*), W1(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv1

call set_symblock_Xaaaav(sleft, W1, nir, nsym, psym) ! -> Xaaaav (allocate) 
call g_sigma_oovv_oovv_no0_x1_type1_eri_v(sv1, iv1, h2_i, Xaaaav, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaaaav)

end subroutine g_if_sigma_oovv_oovv_no0_x1_type1_eri_v



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x1_type1_eri_v(s_v1, i_v1, V2_, W1_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock5), intent(inout) :: W1_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:,:,:)
real*8, allocatable :: Z3_(:,:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o3, i_o3, s_o1, i_o1
integer :: s_o4, i_o4, s_a, i_a
! W1(i,o2,k,o3,v1,a) += (    1.00000000) D3(i,o2,k,o3,o1,o4) V2(v1,a,o1,o4) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(IEOR(s_i,s_o2),s_k) == IEOR(IEOR(s_o3,s_v1),s_a) .and. & 
IEOR(IEOR(s_i,s_o2),s_k) == IEOR(IEOR(s_o3,s_o1),s_o4) .and. &
IEOR(s_v1,s_a) == IEOR(s_o1,s_o4)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o3) > 0 .and. &
   psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o4) > 0) then

! Z1 <-- V2(v1,a,o1,o4) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4)))
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Z1_(i_a, i_o1, i_o4) =  &
  V2_(s_o4, s_o1, s_a)%array(i_o4, i_o1, i_a)
end do
end do
end do
! Z2 <-- D3(i,o2,k,o3,o1,o4) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3)))
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
Z2_(i_o1, i_o4, i_i, i_o2, i_k, i_o3) =  &
  D3_(s_o4, s_o1, s_o3, s_k, s_o2, s_i)%array(i_o4, i_o1, i_o3, i_k, i_o2, i_i)
end do
end do
end do
end do
end do
end do

! Z3 <-- W1(i,o2,k,o3,v1,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o3),&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o4),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o4),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W1(i,o2,k,o3,v1,a)  <-- Z3
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W1_(s_a, s_o3, s_k, s_o2, s_i)%array(i_a, i_o3, i_k, i_o2, i_i) = &
    W1_(s_a, s_o3, s_k, s_o2, s_i)%array(i_a, i_o3, i_k, i_o2, i_i) &
  + Z3_(i_a, i_i, i_o2, i_k, i_o3)
end do
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o3) * &
                psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o4) * 2.0d+00

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

end subroutine g_sigma_oovv_oovv_no0_x1_type1_eri_v



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x1_type1_eri_v(sc, ic, sv1, iv1, T2, W1, S2, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), W1(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv1

call set_symblock_Xaaaav(sleft, W1, nir, nsym, psym) ! -> Xaaaav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x1_type1_eri_v(sc, ic, sv1, iv1, av2_i, Xaaaav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaaaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x1_type1_eri_v



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x1_type1_eri_v(s_c, i_c, s_v1, i_v1, T2_, W1_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
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
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o2, i_o2, s_i, i_i, s_k, i_k, s_a, i_a
! S2(i,k,a,c) += (    2.00000000) T2(o3,o2,c,v1) W1(i,o2,k,o3,v1,a) 
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o3,s_o2) == IEOR(s_c,s_v1) .and. &
IEOR(IEOR(s_i,s_o2),s_k) == IEOR(IEOR(s_o3,s_v1),s_a)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_o3) > 0) then

! Z1 <-- W1(i,o2,k,o3,v1,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
Z1_(i_i, i_k, i_a, i_o2, i_o3) =  &
  W1_(s_a, s_o3, s_k, s_o2, s_i)%array(i_a, i_o3, i_k, i_o2, i_i)
end do
end do
end do
end do
end do
! Z2 <-- T2(o3,o2,c,v1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3)))
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
Z2_(i_o2, i_o3) =  &
  T2_(s_c, s_o2, s_o3)%array(i_c, i_o2, i_o3)
end do
end do

! Z3 <-- S2(i,k,a,c) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a),&
                     1,&
                     psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_o3),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_o3),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a))

! S2(i,k,a,c)  <-- Z3
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) = &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + Z3_(i_i, i_k, i_a)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a) * &
                1 * &
                psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_o3) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no1_x1_type1_eri_v



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x2_type1_eri_v(sc, ic, sv1, iv1, V2, W2, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: V2(*), W2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = IEOR(sc,sv1)

call set_symblock_Xaaaa(sleft, W2, nir, nsym, psym) ! -> Xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x2_type1_eri_v(sc, ic, sv1, iv1, h2_i, Xaaaa, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x2_type1_eri_v



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x2_type1_eri_v(s_c, i_c, s_v1, i_v1, V2_, W2_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock4), intent(inout) :: W2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o3, i_o3, s_k, i_k, s_o2, i_o2, s_o4, i_o4
integer :: s_o1, i_o1
! W2(i,o3,k,o1,c,v1) += (    1.00000000) D3(i,o3,k,o2,o4,o1) V2(c,o2,o4,v1) 
do s_i = 0, nir-1
do s_o3 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o1,s_c),s_v1) .and. & 
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o2,s_o4),s_o1) .and. &
IEOR(s_c,s_o2) == IEOR(s_o4,s_v1)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o1) > 0 .and. &
   psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_o4) > 0) then

! Z1 <-- D3(i,o3,k,o2,o4,o1) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Z1_(i_i, i_o3, i_k, i_o1, i_o2, i_o4) =  &
  D3_(s_o1, s_o4, s_o2, s_k, s_o3, s_i)%array(i_o1, i_o4, i_o2, i_k, i_o3, i_i)
end do
end do
end do
end do
end do
end do
! Z2 <-- V2(c,o2,o4,v1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4)))
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Z2_(i_o2, i_o4) =  &
  V2_(s_v1, s_o4, s_o2)%array(i_v1, i_o4, i_o2)
end do
end do

! Z3 <-- W2(i,o3,k,o1,c,v1) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o1),&
                     1,&
                     psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_o4),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o1),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_o4),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o1))

! W2(i,o3,k,o1,c,v1)  <-- Z3
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W2_(s_o1, s_k, s_o3, s_i)%array(i_o1, i_k, i_o3, i_i) = &
    W2_(s_o1, s_k, s_o3, s_i)%array(i_o1, i_k, i_o3, i_i) &
  + Z3_(i_i, i_o3, i_k, i_o1)
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o1) * &
                1 * &
                psym(I_LENGTH,I_O, s_o2)*psym(I_LENGTH,I_O, s_o4) * 2.0d+00

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

end subroutine g_sigma_oovv_oovv_no0_x2_type1_eri_v



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x2_type1_eri_v(sc, ic, sv1, iv1, T2, W2, S2, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), W2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sc,sv1)

call set_symblock_Xaaaa(sleft, W2, nir, nsym, psym) ! -> Xaaaa (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x2_type1_eri_v(sc, ic, sv1, iv1, av2_i, Xaaaa, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x2_type1_eri_v



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x2_type1_eri_v(s_c, i_c, s_v1, i_v1, T2_, W2_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: W2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o1, i_o1, s_a, i_a, s_i, i_i, s_k, i_k
! S2(i,k,a,c) += (    2.00000000) T2(o3,o1,a,v1) W2(i,o3,k,o1,c,v1) 
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_a = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o3,s_o1) == IEOR(s_a,s_v1) .and. &
IEOR(IEOR(s_i,s_o3),s_k) == IEOR(IEOR(s_o1,s_c),s_v1)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o1) > 0) then

! Z1 <-- T2(o3,o1,a,v1) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1)))
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Z1_(i_a, i_o3, i_o1) =  &
  T2_(s_a, s_o1, s_o3)%array(i_a, i_o1, i_o3)
end do
end do
end do
! Z2 <-- W2(i,o3,k,o1,c,v1) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_o3, i_o1, i_i, i_k) =  &
  W2_(s_o1, s_k, s_o3, s_i)%array(i_o1, i_k, i_o3, i_i)
end do
end do
end do
end do

! Z3 <-- S2(i,k,a,c) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o1),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,k,a,c)  <-- Z3
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) = &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + Z3_(i_a, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no1_x2_type1_eri_v



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x3_type1_eri_v(sv1, iv1, V2, W3, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv1, iv1
real(kind=8), intent(inout) :: V2(*), W3(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sv1

call set_symblock_Xaaaav(sleft, W3, nir, nsym, psym) ! -> Xaaaav (allocate) 
call g_sigma_oovv_oovv_no0_x3_type1_eri_v(sv1, iv1, h2_i, Xaaaav, d3, nir, nsym, psym, flops)

deallocate(h2_i)
deallocate(Xaaaav)

end subroutine g_if_sigma_oovv_oovv_no0_x3_type1_eri_v



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x3_type1_eri_v(s_v1, i_v1, V2_, W3_, D3_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock5), intent(inout) :: W3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:,:,:)
real*8, allocatable :: Z3_(:,:,:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_k, i_k, s_o3, i_o3, s_o1, i_o1
integer :: s_o4, i_o4, s_a, i_a
! W3(i,k,o3,o4,v1,a) += (    1.00000000) D3(i,o2,k,o3,o1,o4) V2(v1,o1,o2,a) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(IEOR(s_i,s_k),s_o3) == IEOR(IEOR(s_o4,s_v1),s_a) .and. & 
IEOR(IEOR(s_i,s_o2),s_k) == IEOR(IEOR(s_o3,s_o1),s_o4) .and. &
IEOR(s_v1,s_o1) == IEOR(s_o2,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4) > 0 .and. &
   psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2) > 0) then

! Z1 <-- V2(v1,o1,o2,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2)))
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Z1_(i_a, i_o1, i_o2) =  &
  V2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
end do
end do
end do
! Z2 <-- D3(i,o2,k,o3,o1,o4) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4)))
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Z2_(i_o1, i_o2, i_i, i_k, i_o3, i_o4) =  &
  D3_(s_o4, s_o1, s_o3, s_k, s_o2, s_i)%array(i_o4, i_o1, i_o3, i_k, i_o2, i_i)
end do
end do
end do
end do
end do
end do

! Z3 <-- W3(i,k,o3,o4,v1,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4),&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W3(i,k,o3,o4,v1,a)  <-- Z3
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
W3_(s_a, s_o4, s_o3, s_k, s_i)%array(i_a, i_o4, i_o3, i_k, i_i) = &
    W3_(s_a, s_o4, s_o3, s_k, s_i)%array(i_a, i_o4, i_o3, i_k, i_i) &
  + Z3_(i_a, i_i, i_k, i_o3, i_o4)
end do
end do
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4) * &
                psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2) * 2.0d+00

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

end subroutine g_sigma_oovv_oovv_no0_x3_type1_eri_v



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x3_type1_eri_v(sc, ic, sv1, iv1, T2, W3, S2, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv1, iv1
real(kind=8), intent(inout) :: T2(*), W3(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sv1

call set_symblock_Xaaaav(sleft, W3, nir, nsym, psym) ! -> Xaaaav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x3_type1_eri_v(sc, ic, sv1, iv1, av2_i, Xaaaav, av2_i2, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(Xaaaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x3_type1_eri_v



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x3_type1_eri_v(s_c, i_c, s_v1, i_v1, T2_, W3_, S2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v1, s_v1
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: W3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:,:,:)
real*8, allocatable :: Z2_(:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_o4, i_o4, s_o3, i_o3, s_i, i_i, s_k, i_k, s_a, i_a
! S2(i,k,a,c) += (    2.00000000) T2(o4,o3,v1,c) W3(i,k,o3,o4,v1,a) 
do s_o4 = 0, nir-1
do s_o3 = 0, nir-1
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_o4,s_o3) == IEOR(s_v1,s_c) .and. &
IEOR(IEOR(s_i,s_k),s_o3) == IEOR(IEOR(s_o4,s_v1),s_a)) then

if(psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4) > 0) then

! Z1 <-- W3(i,k,o3,o4,v1,a) 
allocate(Z1_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4)))
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Z1_(i_i, i_k, i_a, i_o3, i_o4) =  &
  W3_(s_a, s_o4, s_o3, s_k, s_i)%array(i_a, i_o4, i_o3, i_k, i_i)
end do
end do
end do
end do
end do
! Z2 <-- T2(o4,o3,v1,c) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o3):psym(I_END,I_O, s_o3), &
             psym(I_BEGIN,I_O, s_o4):psym(I_END,I_O, s_o4)))
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
Z2_(i_o3, i_o4) =  &
  T2_(s_v1, s_o3, s_o4)%array(i_v1, i_o3, i_o4)
end do
end do

! Z3 <-- S2(i,k,a,c) 
allocate(Z3_(psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k), &
             psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a),&
                     1,&
                     psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a))

! S2(i,k,a,c)  <-- Z3
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) = &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + Z3_(i_i, i_k, i_a)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k)*psym(I_LENGTH,I_V, s_a) * &
                1 * &
                psym(I_LENGTH,I_O, s_o3)*psym(I_LENGTH,I_O, s_o4) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no1_x3_type1_eri_v



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x4_type1_eri_v(sc, ic, sv2, iv2, T2, V2, W4, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sv2, iv2
real(kind=8), intent(inout) :: T2(*), V2(*), W4(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv2, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc

call set_symblock_Xaav(sleft, W4, nir, nsym, psym) ! -> Xaav (allocate) 
call g_sigma_oovv_oovv_no0_x4_type1_eri_v(sc, ic, sv2, iv2, av2_i, h2_i, Xaav, nir, nsym, psym, flops)

deallocate(av2_i)
deallocate(h2_i)
deallocate(Xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x4_type1_eri_v



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no0_x4_type1_eri_v(s_c, i_c, s_v2, i_v2, T2_, V2_, W4_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_v2, s_v2
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: W4_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:)
real*8, allocatable :: Z2_(:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1, s_a, i_a
! W4(o1,o2,c,a) += (    1.00000000) T2(o1,o2,v1,v2) V2(c,v2,a,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_c,s_a) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_v1,s_v2) .and. &
IEOR(s_c,s_v2) == IEOR(s_a,s_v1)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2) > 0 .and. &
   psym(I_LENGTH,I_V, s_v1) > 0) then

! Z1 <-- V2(c,v2,a,v1) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_V, s_v1):psym(I_END,I_V, s_v1)))
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Z1_(i_a, i_v1) =  &
  V2_(s_v1, s_a, s_v2)%array(i_v1, i_a, i_v2)
end do
end do
! Z2 <-- T2(o1,o2,v1,v2) 
allocate(Z2_(psym(I_BEGIN,I_V, s_v1):psym(I_END,I_V, s_v1), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2)))
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Z2_(i_v1, i_o1, i_o2) =  &
  T2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1)
end do
end do
end do

! Z3 <-- W4(o1,o2,c,a) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2),&
                     psym(I_LENGTH,I_V, s_v1),&
                     1.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_V, s_v1),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! W4(o1,o2,c,a)  <-- Z3
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
W4_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1) = &
    W4_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1) &
  + Z3_(i_a, i_o1, i_o2)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2) * &
                psym(I_LENGTH,I_V, s_v1) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no0_x4_type1_eri_v



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x4_type1_eri_v(sc, ic, W4, S2, nir, nsym, psym, flops)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: W4(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops 
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = sc

call set_symblock_Xaav(sleft, W4, nir, nsym, psym) ! -> Xaav (allocate) 
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x4_type1_eri_v(sc, ic, Xaav, av2_i2, d2, nir, nsym, psym, flops)

deallocate(Xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x4_type1_eri_v



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_oovv_no1_x4_type1_eri_v(s_c, i_c, W4_, S2_, D2_, nir, nsym, psym, flops)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6

implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Flop count
real*8, intent(inout) :: flops

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: W4_(0:nir-1, 0:nir-1, 0:nir-1)

! Intermediate arrays                
real*8, allocatable :: Z1_(:,:,:)
real*8, allocatable :: Z2_(:,:,:,:)
real*8, allocatable :: Z3_(:,:,:)
! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2, s_a, i_a
! S2(i,k,a,c) += (    2.00000000) D2(i,o1,k,o2) W4(o1,o2,c,a) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_o2) .and. &
IEOR(s_o1,s_o2) == IEOR(s_c,s_a)) then

if(psym(I_LENGTH,I_V, s_a) > 0 .and. &
   psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) > 0 .and. &
   psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2) > 0) then

! Z1 <-- W4(o1,o2,c,a) 
allocate(Z1_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2)))
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Z1_(i_a, i_o1, i_o2) =  &
  W4_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
end do
end do
end do
! Z2 <-- D2(i,o1,k,o2) 
allocate(Z2_(psym(I_BEGIN,I_O, s_o1):psym(I_END,I_O, s_o1), &
             psym(I_BEGIN,I_O, s_o2):psym(I_END,I_O, s_o2), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Z2_(i_o1, i_o2, i_i, i_k) =  &
  D2_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i)
end do
end do
end do
end do

! Z3 <-- S2(i,k,a,c) 
allocate(Z3_(psym(I_BEGIN,I_V, s_a):psym(I_END,I_V, s_a), &
             psym(I_BEGIN,I_O, s_i):psym(I_END,I_O, s_i), &
             psym(I_BEGIN,I_O, s_k):psym(I_END,I_O, s_k)))

! Gemm Z1 * Z2 to form Z3
call dgemm('n', 'n', psym(I_LENGTH,I_V, s_a),&
                     psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k),&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2),&
                     2.00000000d+00, &
                     Z1_,&
                     psym(I_LENGTH,I_V, s_a),&
                     Z2_,&
                     psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2),&
                     0.0d+00,&
                     Z3_,&
                     psym(I_LENGTH,I_V, s_a))

! S2(i,k,a,c)  <-- Z3
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) = &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + Z3_(i_a, i_i, i_k)
end do
end do
end do

! Flop count
flops = flops + psym(I_LENGTH,I_V, s_a) * &
                psym(I_LENGTH,I_O, s_i)*psym(I_LENGTH,I_O, s_k) * &
                psym(I_LENGTH,I_O, s_o1)*psym(I_LENGTH,I_O, s_o2) * 2.0d+00

deallocate(Z1_, Z2_, Z3_)
end if ! Dim Const

end if ! Irrep Cond
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_oovv_oovv_no1_x4_type1_eri_v

