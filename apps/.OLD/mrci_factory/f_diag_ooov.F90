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
subroutine g_if_diag_ooov_no0_x0(sa, ia, h, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: h(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x0(sa, ia, h1, av2_i2, d3, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x0(s_a, i_a, h_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D3(i,m,k,k,o1,i) h(m,o1) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_k,s_o1),s_i) .and. &
IEOR(s_m,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_i, s_o1, s_k, s_k, s_m, s_i)%array(i_i, i_o1, i_k, i_k, i_m, i_i) & 
  * h_(s_o1, s_m)%array(i_o1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x1(sa, ia, h, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: h(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x1(sa, ia, h1, av2_i2, d2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x1(s_a, i_a, h_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D2(i,i,k,k) h(m,m) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k) .and. &
IEOR(s_m,s_m) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i) & 
  * h_(s_m, s_m)%array(i_m, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x2(sa, ia, h, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: h(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x2(sa, ia, h1, av2_i2, d3, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x2(s_a, i_a, h_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D3(i,m,k,k,m,o1) h(i,o1) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_k,s_m),s_o1) .and. &
IEOR(s_i,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_o1, s_m, s_k, s_k, s_m, s_i)%array(i_o1, i_m, i_k, i_k, i_m, i_i) & 
  * h_(s_o1, s_i)%array(i_o1, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x3(sa, ia, h, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: h(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x3(sa, ia, h1, av2_i2, d2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x3(s_a, i_a, h_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D2(i,o1,k,k) h(i,o1) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_k) .and. &
IEOR(s_i,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_k, s_k, s_o1, s_i)%array(i_k, i_k, i_o1, i_i) & 
  * h_(s_o1, s_i)%array(i_o1, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x4(sa, ia, h, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: h(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x4(sa, ia, h1, av2_i2, d3, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x4(s_a, i_a, h_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D3(i,m,k,o1,m,i) h(k,o1) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o1,s_m),s_i) .and. &
IEOR(s_k,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_i, s_m, s_o1, s_k, s_m, s_i)%array(i_i, i_m, i_o1, i_k, i_m, i_i) & 
  * h_(s_o1, s_k)%array(i_o1, i_k)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x5(sa, ia, h, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: h(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x5(sa, ia, h1, av2_i2, d2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x5



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x5(s_a, i_a, h_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D2(i,i,k,o1) h(k,o1) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_o1) .and. &
IEOR(s_k,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_i, s_i)%array(i_o1, i_k, i_i, i_i) & 
  * h_(s_o1, s_k)%array(i_o1, i_k)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x5



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x6(sa, ia, h, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: h(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x6(sa, ia, h1, av2_i2, d3, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x6



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x6(s_a, i_a, h_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D3(i,m,k,k,m,i) h(a,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_k,s_m),s_i) .and. &
IEOR(s_a,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_i, s_m, s_k, s_k, s_m, s_i)%array(i_i, i_m, i_k, i_k, i_m, i_i) & 
  * h_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x6



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x7(sa, ia, h, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: h(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x7(sa, ia, h1, av2_i2, d2, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x7



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x7(s_a, i_a, h_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D2(i,i,k,k) h(a,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k) .and. &
IEOR(s_a,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i) & 
  * h_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x7



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y8(sc1, ic1, V2, Y0, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y0(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y0, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ooov_y8(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ooov_y8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y8(s_c1, i_c1, V2_, Y0_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y0_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_m, i_m, s_o1, i_o1
! Y0(m,o1)  <-- 
! (    1.00000000)  V2(c1,c1,m,o1) 
do s_m = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_m,s_o1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_m,s_o1)) then
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Y0_(s_o1, s_m)%array(i_o1, i_m) =  &
    Y0_(s_o1, s_m)%array(i_o1, i_m) &
  + 1.00000000d+00 & 
  * V2_(s_o1, s_m, s_c1)%array(i_o1, i_m, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x8(sa, ia, Y0, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y0(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y0, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x8(sa, ia, Yaa, av2_i2, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x8



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x8(s_a, i_a, Y0_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y0_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! Hdiag(i,k,m,a)  <-- 
! (    2.00000000) D3(i,m,k,k,o1,i) Y0(m,o1) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_k,s_o1),s_i) .and. &
IEOR(s_m,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 2.00000000d+00 & 
  * D3_(s_i, s_o1, s_k, s_k, s_m, s_i)%array(i_i, i_o1, i_k, i_k, i_m, i_i) & 
  * Y0_(s_o1, s_m)%array(i_o1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x8



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y9(sc1, ic1, V2, Y1, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y1(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y1, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ooov_y9(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ooov_y9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y9(s_c1, i_c1, V2_, Y1_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y1_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_m, i_m, s_o1, i_o1
! Y1(m,o1)  <-- 
! (    1.00000000)  V2(c1,m,c1,o1) 
do s_m = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_m,s_o1) == 0 .and. & 
IEOR(s_c1,s_m) == IEOR(s_c1,s_o1)) then
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Y1_(s_o1, s_m)%array(i_o1, i_m) =  &
    Y1_(s_o1, s_m)%array(i_o1, i_m) &
  + 1.00000000d+00 & 
  * V2_(s_o1, s_c1, s_m)%array(i_o1, i_c1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x9(sa, ia, Y1, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y1(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y1, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x9(sa, ia, Yaa, av2_i2, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x9



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x9(s_a, i_a, Y1_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y1_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D3(i,m,k,k,o1,i) Y1(m,o1) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_k,s_o1),s_i) .and. &
IEOR(s_m,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_i, s_o1, s_k, s_k, s_m, s_i)%array(i_i, i_o1, i_k, i_k, i_m, i_i) & 
  * Y1_(s_o1, s_m)%array(i_o1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x9



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x10(sa, ia, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x10(sa, ia, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x10



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x10(s_a, i_a, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! Hdiag(i,k,m,a)  <-- 
! (    2.00000000) D3(i,m,k,k,o1,i) V2(a,a,m,o1) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_k,s_o1),s_i) .and. &
IEOR(s_a,s_a) == IEOR(s_m,s_o1)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 2.00000000d+00 & 
  * D3_(s_i, s_o1, s_k, s_k, s_m, s_i)%array(i_i, i_o1, i_k, i_k, i_m, i_i) & 
  * V2_(s_o1, s_m, s_a)%array(i_o1, i_m, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x10



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y11(sc1, ic1, V2, Y2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y2, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ooov_y11(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ooov_y11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y11(s_c1, i_c1, V2_, Y2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y2_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_m, i_m
! Y2(m,m)  <-- 
! (    1.00000000)  V2(c1,c1,m,m) 
do s_m = 0, nir-1
if( &
IEOR(s_m,s_m) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_m,s_m)) then
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Y2_(s_m, s_m)%array(i_m, i_m) =  &
    Y2_(s_m, s_m)%array(i_m, i_m) &
  + 1.00000000d+00 & 
  * V2_(s_m, s_m, s_c1)%array(i_m, i_m, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x11(sa, ia, Y2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y2, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x11(sa, ia, Yaa, av2_i2, d2, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x11



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x11(s_a, i_a, Y2_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y2_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (    2.00000000) D2(i,i,k,k) Y2(m,m) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k) .and. &
IEOR(s_m,s_m) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 2.00000000d+00 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i) & 
  * Y2_(s_m, s_m)%array(i_m, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x11



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y12(sc1, ic1, V2, Y3, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y3(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y3, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ooov_y12(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ooov_y12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y12(s_c1, i_c1, V2_, Y3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y3_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_m, i_m
! Y3(m,m)  <-- 
! (    1.00000000)  V2(c1,m,c1,m) 
do s_m = 0, nir-1
if( &
IEOR(s_m,s_m) == 0 .and. & 
IEOR(s_c1,s_m) == IEOR(s_c1,s_m)) then
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Y3_(s_m, s_m)%array(i_m, i_m) =  &
    Y3_(s_m, s_m)%array(i_m, i_m) &
  + 1.00000000d+00 & 
  * V2_(s_m, s_c1, s_m)%array(i_m, i_c1, i_m)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x12(sa, ia, Y3, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y3(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y3, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x12(sa, ia, Yaa, av2_i2, d2, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x12



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x12(s_a, i_a, Y3_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y3_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D2(i,i,k,k) Y3(m,m) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k) .and. &
IEOR(s_m,s_m) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i) & 
  * Y3_(s_m, s_m)%array(i_m, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x12



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y13(sc1, ic1, V2, Y4, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y4(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y4, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ooov_y13(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ooov_y13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y13(s_c1, i_c1, V2_, Y4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y4_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a
! Y4(a,a)  <-- 
! (    1.00000000)  V2(c1,c1,a,a) 
do s_a = 0, nir-1
if( &
IEOR(s_a,s_a) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_a)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y4_(s_a, s_a)%array(i_a, i_a) =  &
    Y4_(s_a, s_a)%array(i_a, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_a, s_c1)%array(i_a, i_a, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x13(sa, ia, Y4, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y4(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y4, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x13(sa, ia, Yvv, av2_i2, d2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x13



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x13(s_a, i_a, Y4_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y4_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (    2.00000000) D2(i,i,k,k) Y4(a,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k) .and. &
IEOR(s_a,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 2.00000000d+00 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i) & 
  * Y4_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x13



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y14(sc1, ic1, V2, Y5, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y5(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y5, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ooov_y14(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ooov_y14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y14(s_c1, i_c1, V2_, Y5_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y5_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a
! Y5(a,a)  <-- 
! (    1.00000000)  V2(c1,a,c1,a) 
do s_a = 0, nir-1
if( &
IEOR(s_a,s_a) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_a)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y5_(s_a, s_a)%array(i_a, i_a) =  &
    Y5_(s_a, s_a)%array(i_a, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_c1, s_a)%array(i_a, i_c1, i_a)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x14(sa, ia, Y5, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y5(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y5, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x14(sa, ia, Yvv, av2_i2, d2, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x14



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x14(s_a, i_a, Y5_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y5_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D2(i,i,k,k) Y5(a,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k) .and. &
IEOR(s_a,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i) & 
  * Y5_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x14



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x15(sa, ia, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x15(sa, ia, h2_i, av2_i2, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x15



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x15(s_a, i_a, V2_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D2(i,i,k,k) V2(a,a,m,m) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k) .and. &
IEOR(s_a,s_a) == IEOR(s_m,s_m)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i) & 
  * V2_(s_m, s_m, s_a)%array(i_m, i_m, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x15



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y16(sc1, ic1, V2, Y6, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y6(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y6, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ooov_y16(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ooov_y16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y16(s_c1, i_c1, V2_, Y6_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y6_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1
! Y6(i,o1)  <-- 
! (    1.00000000)  V2(c1,c1,i,o1) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_o1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_i,s_o1)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Y6_(s_o1, s_i)%array(i_o1, i_i) =  &
    Y6_(s_o1, s_i)%array(i_o1, i_i) &
  + 1.00000000d+00 & 
  * V2_(s_o1, s_i, s_c1)%array(i_o1, i_i, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x16(sa, ia, Y6, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y6(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y6, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x16(sa, ia, Yaa, av2_i2, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x16



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x16(s_a, i_a, Y6_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y6_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! Hdiag(i,k,m,a)  <-- 
! (   -2.00000000) D3(i,m,k,k,m,o1) Y6(i,o1) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_k,s_m),s_o1) .and. &
IEOR(s_i,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 2.00000000d+00 & 
  * D3_(s_o1, s_m, s_k, s_k, s_m, s_i)%array(i_o1, i_m, i_k, i_k, i_m, i_i) & 
  * Y6_(s_o1, s_i)%array(i_o1, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x16



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y17(sc1, ic1, V2, Y7, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y7(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y7, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ooov_y17(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ooov_y17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y17(s_c1, i_c1, V2_, Y7_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y7_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1
! Y7(i,o1)  <-- 
! (    1.00000000)  V2(c1,i,c1,o1) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_o1) == 0 .and. & 
IEOR(s_c1,s_i) == IEOR(s_c1,s_o1)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Y7_(s_o1, s_i)%array(i_o1, i_i) =  &
    Y7_(s_o1, s_i)%array(i_o1, i_i) &
  + 1.00000000d+00 & 
  * V2_(s_o1, s_c1, s_i)%array(i_o1, i_c1, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x17(sa, ia, Y7, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y7(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y7, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x17(sa, ia, Yaa, av2_i2, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x17



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x17(s_a, i_a, Y7_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y7_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D3(i,m,k,k,m,o1) Y7(i,o1) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_k,s_m),s_o1) .and. &
IEOR(s_i,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o1, s_m, s_k, s_k, s_m, s_i)%array(i_o1, i_m, i_k, i_k, i_m, i_i) & 
  * Y7_(s_o1, s_i)%array(i_o1, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x17



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y18(sc1, ic1, V2, Y8, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y8(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y8, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ooov_y18(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ooov_y18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y18(s_c1, i_c1, V2_, Y8_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y8_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1
! Y8(i,o1)  <-- 
! (    1.00000000)  V2(c1,c1,i,o1) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_o1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_i,s_o1)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Y8_(s_o1, s_i)%array(i_o1, i_i) =  &
    Y8_(s_o1, s_i)%array(i_o1, i_i) &
  + 1.00000000d+00 & 
  * V2_(s_o1, s_i, s_c1)%array(i_o1, i_i, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x18(sa, ia, Y8, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y8(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y8, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x18(sa, ia, Yaa, av2_i2, d2, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x18



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x18(s_a, i_a, Y8_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y8_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (   -2.00000000) D2(i,o1,k,k) Y8(i,o1) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_k) .and. &
IEOR(s_i,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 2.00000000d+00 & 
  * D2_(s_k, s_k, s_o1, s_i)%array(i_k, i_k, i_o1, i_i) & 
  * Y8_(s_o1, s_i)%array(i_o1, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x18



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y19(sc1, ic1, V2, Y9, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y9(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y9, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ooov_y19(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ooov_y19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y19(s_c1, i_c1, V2_, Y9_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y9_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1
! Y9(i,o1)  <-- 
! (    1.00000000)  V2(c1,i,c1,o1) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_o1) == 0 .and. & 
IEOR(s_c1,s_i) == IEOR(s_c1,s_o1)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Y9_(s_o1, s_i)%array(i_o1, i_i) =  &
    Y9_(s_o1, s_i)%array(i_o1, i_i) &
  + 1.00000000d+00 & 
  * V2_(s_o1, s_c1, s_i)%array(i_o1, i_c1, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x19(sa, ia, Y9, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y9(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y9, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x19(sa, ia, Yaa, av2_i2, d2, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x19



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x19(s_a, i_a, Y9_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y9_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D2(i,o1,k,k) Y9(i,o1) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_k) .and. &
IEOR(s_i,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_k, s_k, s_o1, s_i)%array(i_k, i_k, i_o1, i_i) & 
  * Y9_(s_o1, s_i)%array(i_o1, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x19



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y20(sc1, ic1, V2, Y10, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y10(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y10, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ooov_y20(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ooov_y20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y20(s_c1, i_c1, V2_, Y10_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y10_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_o1, i_o1
! Y10(k,o1)  <-- 
! (    1.00000000)  V2(c1,c1,k,o1) 
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_k,s_o1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_k,s_o1)) then
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Y10_(s_o1, s_k)%array(i_o1, i_k) =  &
    Y10_(s_o1, s_k)%array(i_o1, i_k) &
  + 1.00000000d+00 & 
  * V2_(s_o1, s_k, s_c1)%array(i_o1, i_k, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x20(sa, ia, Y10, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y10(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y10, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x20(sa, ia, Yaa, av2_i2, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x20



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x20(s_a, i_a, Y10_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y10_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! Hdiag(i,k,m,a)  <-- 
! (   -2.00000000) D3(i,m,k,o1,m,i) Y10(k,o1) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o1,s_m),s_i) .and. &
IEOR(s_k,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 2.00000000d+00 & 
  * D3_(s_i, s_m, s_o1, s_k, s_m, s_i)%array(i_i, i_m, i_o1, i_k, i_m, i_i) & 
  * Y10_(s_o1, s_k)%array(i_o1, i_k)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x20



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y21(sc1, ic1, V2, Y11, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y11(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y11, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ooov_y21(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ooov_y21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y21(s_c1, i_c1, V2_, Y11_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y11_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_o1, i_o1
! Y11(k,o1)  <-- 
! (    1.00000000)  V2(c1,k,c1,o1) 
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_k,s_o1) == 0 .and. & 
IEOR(s_c1,s_k) == IEOR(s_c1,s_o1)) then
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Y11_(s_o1, s_k)%array(i_o1, i_k) =  &
    Y11_(s_o1, s_k)%array(i_o1, i_k) &
  + 1.00000000d+00 & 
  * V2_(s_o1, s_c1, s_k)%array(i_o1, i_c1, i_k)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x21(sa, ia, Y11, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y11(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y11, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x21(sa, ia, Yaa, av2_i2, d3, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x21



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x21(s_a, i_a, Y11_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y11_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k, s_o1, i_o1
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D3(i,m,k,o1,m,i) Y11(k,o1) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o1,s_m),s_i) .and. &
IEOR(s_k,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_i, s_m, s_o1, s_k, s_m, s_i)%array(i_i, i_m, i_o1, i_k, i_m, i_i) & 
  * Y11_(s_o1, s_k)%array(i_o1, i_k)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x21



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y22(sc1, ic1, V2, Y12, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y12(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y12, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ooov_y22(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ooov_y22



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y22(s_c1, i_c1, V2_, Y12_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y12_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_o1, i_o1
! Y12(k,o1)  <-- 
! (    1.00000000)  V2(c1,c1,k,o1) 
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_k,s_o1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_k,s_o1)) then
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Y12_(s_o1, s_k)%array(i_o1, i_k) =  &
    Y12_(s_o1, s_k)%array(i_o1, i_k) &
  + 1.00000000d+00 & 
  * V2_(s_o1, s_k, s_c1)%array(i_o1, i_k, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y22



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x22(sa, ia, Y12, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y12(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y12, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x22(sa, ia, Yaa, av2_i2, d2, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x22



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x22(s_a, i_a, Y12_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y12_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (   -2.00000000) D2(i,i,k,o1) Y12(k,o1) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_o1) .and. &
IEOR(s_k,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 2.00000000d+00 & 
  * D2_(s_o1, s_k, s_i, s_i)%array(i_o1, i_k, i_i, i_i) & 
  * Y12_(s_o1, s_k)%array(i_o1, i_k)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x22



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y23(sc1, ic1, V2, Y13, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y13(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yaa(sleft2, Y13, nir, nsym, psym) ! -> Yaa (allocate) 
call g_diag_ooov_y23(sc1, ic1, h2_i, Yaa, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yaa)

end subroutine g_if_diag_ooov_y23



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y23(s_c1, i_c1, V2_, Y13_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y13_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_o1, i_o1
! Y13(k,o1)  <-- 
! (    1.00000000)  V2(c1,k,c1,o1) 
do s_k = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_k,s_o1) == 0 .and. & 
IEOR(s_c1,s_k) == IEOR(s_c1,s_o1)) then
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Y13_(s_o1, s_k)%array(i_o1, i_k) =  &
    Y13_(s_o1, s_k)%array(i_o1, i_k) &
  + 1.00000000d+00 & 
  * V2_(s_o1, s_c1, s_k)%array(i_o1, i_c1, i_k)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y23



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x23(sa, ia, Y13, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y13(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yaa(sleft2, Y13, nir, nsym, psym) ! -> Yaa (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x23(sa, ia, Yaa, av2_i2, d2, nir, nsym, psym)

deallocate(Yaa)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x23



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x23(s_a, i_a, Y13_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y13_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D2(i,i,k,o1) Y13(k,o1) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_o1) .and. &
IEOR(s_k,s_o1) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_i, s_i)%array(i_o1, i_k, i_i, i_i) & 
  * Y13_(s_o1, s_k)%array(i_o1, i_k)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x23



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y24(sc1, ic1, V2, Y14, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y14(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y14, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ooov_y24(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ooov_y24



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y24(s_c1, i_c1, V2_, Y14_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y14_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a
! Y14(a,a)  <-- 
! (    1.00000000)  V2(c1,c1,a,a) 
do s_a = 0, nir-1
if( &
IEOR(s_a,s_a) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_a,s_a)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y14_(s_a, s_a)%array(i_a, i_a) =  &
    Y14_(s_a, s_a)%array(i_a, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_a, s_c1)%array(i_a, i_a, i_c1)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y24



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x24(sa, ia, Y14, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y14(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y14, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x24(sa, ia, Yvv, av2_i2, d3, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x24



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x24(s_a, i_a, Y14_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y14_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k
! Hdiag(i,k,m,a)  <-- 
! (    2.00000000) D3(i,m,k,k,m,i) Y14(a,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_k,s_m),s_i) .and. &
IEOR(s_a,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 2.00000000d+00 & 
  * D3_(s_i, s_m, s_k, s_k, s_m, s_i)%array(i_i, i_m, i_k, i_k, i_m, i_i) & 
  * Y14_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x24



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_y25(sc1, ic1, V2, Y15, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc1, ic1
real(kind=8), intent(inout) :: V2(*), Y15(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft2 = 0
call set_symblock_Yvv(sleft2, Y15, nir, nsym, psym) ! -> Yvv (allocate) 
call g_diag_ooov_y25(sc1, ic1, h2_i, Yvv, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yvv)

end subroutine g_if_diag_ooov_y25



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_y25(s_c1, i_c1, V2_, Y15_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c1, s_c1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: Y15_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a
! Y15(a,a)  <-- 
! (    1.00000000)  V2(c1,a,c1,a) 
do s_a = 0, nir-1
if( &
IEOR(s_a,s_a) == 0 .and. & 
IEOR(s_c1,s_a) == IEOR(s_c1,s_a)) then
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
Y15_(s_a, s_a)%array(i_a, i_a) =  &
    Y15_(s_a, s_a)%array(i_a, i_a) &
  + 1.00000000d+00 & 
  * V2_(s_a, s_c1, s_a)%array(i_a, i_c1, i_a)
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_y25



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x25(sa, ia, Y15, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Y15(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft2 = 0
call set_symblock_Yvv(sleft2, Y15, nir, nsym, psym) ! -> Yvv (allocate) 
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x25(sa, ia, Yvv, av2_i2, d3, nir, nsym, psym)

deallocate(Yvv)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x25



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x25(s_a, i_a, Y15_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y15_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D3(i,m,k,k,m,i) Y15(a,a) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_k,s_m),s_i) .and. &
IEOR(s_a,s_a) == 0) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_i, s_m, s_k, s_k, s_m, s_i)%array(i_i, i_m, i_k, i_k, i_m, i_i) & 
  * Y15_(s_a, s_a)%array(i_a, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x25



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x26(sa, ia, si, ii, sm, im, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, si, ii, sm, im
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sm, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x26(sa, ia, si, ii, sm, im, h2_i, av2_i2, d4_ij, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x26



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x26(s_a, i_a, s_i, i_i, s_m, i_m, V2_, Hdiag_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_m, s_m
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_o1, i_o1, s_o3, i_o3, s_o2, i_o2
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D4(m,i,k,k,o3,o1,i,o2) V2(m,o2,o1,o3) 
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_m,s_i),IEOR(s_k,s_k)) == IEOR(IEOR(s_o3,s_o1),IEOR(s_i,s_o2)) .and. &
IEOR(s_m,s_o2) == IEOR(s_o1,s_o3)) then
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D4_(s_o2, s_i, s_o1, s_o3, s_k, s_k)%array(i_o2, i_i, i_o1, i_o3, i_k, i_k) & 
  * V2_(s_o3, s_o1, s_o2)%array(i_o3, i_o1, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x26



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x27(sa, ia, sm, im, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sm, im
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sm, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x27(sa, ia, sm, im, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x27



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x27(s_a, i_a, s_m, i_m, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_m, s_m
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D3(i,i,k,k,o1,o2) V2(m,m,o1,o2) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_i),s_k) == IEOR(IEOR(s_k,s_o1),s_o2) .and. &
IEOR(s_m,s_m) == IEOR(s_o1,s_o2)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o1, s_k, s_k, s_i, s_i)%array(i_o2, i_o1, i_k, i_k, i_i, i_i) & 
  * V2_(s_o2, s_o1, s_m)%array(i_o2, i_o1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x27



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x28(sa, ia, sm, im, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sm, im
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sm, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x28(sa, ia, sm, im, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x28



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x28(s_a, i_a, s_m, i_m, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_m, s_m
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D3(i,o1,k,k,o2,i) V2(m,o1,m,o2) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_o1),s_k) == IEOR(IEOR(s_k,s_o2),s_i) .and. &
IEOR(s_m,s_o1) == IEOR(s_m,s_o2)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_i, s_o2, s_k, s_k, s_o1, s_i)%array(i_i, i_o2, i_k, i_k, i_o1, i_i) & 
  * V2_(s_o2, s_m, s_o1)%array(i_o2, i_m, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x28



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x29(sa, ia, si, ii, sm, im, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, si, ii, sm, im
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(si, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x29(sa, ia, si, ii, sm, im, h2_i, av2_i2, d4_ij, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x29



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x29(s_a, i_a, s_i, i_i, s_m, i_m, V2_, Hdiag_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_m, s_m
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_o2, i_o2, s_o3, i_o3, s_o1, i_o1
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D4(i,m,k,k,m,o2,o3,o1) V2(i,o2,o1,o3) 
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),IEOR(s_k,s_k)) == IEOR(IEOR(s_m,s_o2),IEOR(s_o3,s_o1)) .and. &
IEOR(s_i,s_o2) == IEOR(s_o1,s_o3)) then
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D4_(s_o1, s_o3, s_o2, s_m, s_k, s_k)%array(i_o1, i_o3, i_o2, i_m, i_k, i_k) & 
  * V2_(s_o3, s_o1, s_o2)%array(i_o3, i_o1, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x29



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x30(sa, ia, si, ii, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, si, ii
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(si, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x30(sa, ia, si, ii, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x30



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x30(s_a, i_a, s_i, i_i, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_k, i_k, s_o3, i_o3, s_o1, i_o1, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D3(i,o2,k,k,o3,o1) V2(i,o2,o1,o3) 
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_o2),s_k) == IEOR(IEOR(s_k,s_o3),s_o1) .and. &
IEOR(s_i,s_o2) == IEOR(s_o1,s_o3)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_o1, s_o3, s_k, s_k, s_o2, s_i)%array(i_o1, i_o3, i_k, i_k, i_o2, i_i) & 
  * V2_(s_o3, s_o1, s_o2)%array(i_o3, i_o1, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x30



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x31(sa, ia, si, ii, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, si, ii
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(si, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x31(sa, ia, si, ii, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x31



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x31(s_a, i_a, s_i, i_i, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_k, i_k, s_m, i_m, s_o1, i_o1
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D3(i,o2,k,k,m,o1) V2(i,o1,m,o2) 
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_m = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_o2),s_k) == IEOR(IEOR(s_k,s_m),s_o1) .and. &
IEOR(s_i,s_o1) == IEOR(s_m,s_o2)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_o1, s_m, s_k, s_k, s_o2, s_i)%array(i_o1, i_m, i_k, i_k, i_o2, i_i) & 
  * V2_(s_o2, s_m, s_o1)%array(i_o2, i_m, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x31



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x32(sa, ia, sk, ik, so2, io2, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sk, ik, so2, io2
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sk, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x32(sa, ia, sk, ik, so2, io2, h2_i, av2_i2, d4_ij, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x32



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x32(s_a, i_a, s_k, i_k, s_o2, i_o2, V2_, Hdiag_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_k, s_k
integer, intent(in) :: i_o2, s_o2
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_o3, i_o3, s_o1, i_o1
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D4(k,o2,i,m,m,i,o3,o1) V2(k,o2,o1,o3) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_k,s_o2),IEOR(s_i,s_m)) == IEOR(IEOR(s_m,s_i),IEOR(s_o3,s_o1)) .and. &
IEOR(s_k,s_o2) == IEOR(s_o1,s_o3)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D4_(s_o1, s_o3, s_i, s_m, s_m, s_i)%array(i_o1, i_o3, i_i, i_m, i_m, i_i) & 
  * V2_(s_o3, s_o1, s_o2)%array(i_o3, i_o1, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x32



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x33(sa, ia, sk, ik, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sk, ik
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sk, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x33(sa, ia, sk, ik, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x33



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x33(s_a, i_a, s_k, i_k, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_k, s_k
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o2, i_o2, s_o3, i_o3, s_o1, i_o1, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D3(i,i,k,o2,o3,o1) V2(k,o2,o1,o3) 
do s_i = 0, nir-1
do s_o2 = 0, nir-1
do s_o3 = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_i),s_k) == IEOR(IEOR(s_o2,s_o3),s_o1) .and. &
IEOR(s_k,s_o2) == IEOR(s_o1,s_o3)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_o1, s_o3, s_o2, s_k, s_i, s_i)%array(i_o1, i_o3, i_o2, i_k, i_i, i_i) & 
  * V2_(s_o3, s_o1, s_o2)%array(i_o3, i_o1, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x33



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x34(sa, ia, sk, ik, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sk, ik
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sk, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x34(sa, ia, sk, ik, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x34



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x34(s_a, i_a, s_k, i_k, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_k, s_k
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_o1, i_o1, s_o2, i_o2
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D3(i,m,o1,k,o2,i) V2(k,o1,m,o2) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_o1) == IEOR(IEOR(s_k,s_o2),s_i) .and. &
IEOR(s_k,s_o1) == IEOR(s_m,s_o2)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_i, s_o2, s_k, s_o1, s_m, s_i)%array(i_i, i_o2, i_k, i_o1, i_m, i_i) & 
  * V2_(s_o2, s_m, s_o1)%array(i_o2, i_m, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x34



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x35(sa, ia, si, ii, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, si, ii
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(si, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x35(sa, ia, si, ii, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x35



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x35(s_a, i_a, s_i, i_i, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_m, i_m, s_k, i_k, s_o2, i_o2, s_o1, i_o1
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D3(i,m,k,o2,m,o1) V2(i,o1,k,o2) 
do s_m = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_o2,s_m),s_o1) .and. &
IEOR(s_i,s_o1) == IEOR(s_k,s_o2)) then
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D3_(s_o1, s_m, s_o2, s_k, s_m, s_i)%array(i_o1, i_m, i_o2, i_k, i_m, i_i) & 
  * V2_(s_o2, s_k, s_o1)%array(i_o2, i_k, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x35



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x36(sa, ia, si, ii, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, si, ii
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(si, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x36(sa, ia, si, ii, h2_i, av2_i2, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x36



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x36(s_a, i_a, s_i, i_i, V2_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_k, i_k, s_o1, i_o1, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (   -1.00000000) D2(i,o2,k,o1) V2(i,o2,k,o1) 
do s_o2 = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_o2) == IEOR(s_k,s_o1) .and. &
IEOR(s_i,s_o2) == IEOR(s_k,s_o1)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  - 1.00000000d+00 & 
  * D2_(s_o1, s_k, s_o2, s_i)%array(i_o1, i_k, i_o2, i_i) & 
  * V2_(s_o1, s_k, s_o2)%array(i_o1, i_k, i_o2)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x36



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x37(sa, ia, si, ii, sm, im, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, si, ii, sm, im
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x37(sa, ia, si, ii, sm, im, h2_i, av2_i2, d4_ij, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x37



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x37(s_a, i_a, s_i, i_i, s_m, i_m, V2_, Hdiag_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_m, s_m
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_o1, i_o1, s_o2, i_o2
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D4(i,m,k,k,m,i,o1,o2) V2(a,a,o1,o2) 
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),IEOR(s_k,s_k)) == IEOR(IEOR(s_m,s_i),IEOR(s_o1,s_o2)) .and. &
IEOR(s_a,s_a) == IEOR(s_o1,s_o2)) then
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D4_(s_o2, s_o1, s_i, s_m, s_k, s_k)%array(i_o2, i_o1, i_i, i_m, i_k, i_k) & 
  * V2_(s_o2, s_o1, s_a)%array(i_o2, i_o1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x37



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x38(sa, ia, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x38(sa, ia, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x38



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x38(s_a, i_a, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D3(i,i,k,k,o1,o2) V2(a,a,o1,o2) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_i),s_k) == IEOR(IEOR(s_k,s_o1),s_o2) .and. &
IEOR(s_a,s_a) == IEOR(s_o1,s_o2)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_o2, s_o1, s_k, s_k, s_i, s_i)%array(i_o2, i_o1, i_k, i_k, i_i, i_i) & 
  * V2_(s_o2, s_o1, s_a)%array(i_o2, i_o1, i_a)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x38



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x39(sa, ia, si, ii, sm, im, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, si, ii, sm, im
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x39(sa, ia, si, ii, sm, im, h2_i, av2_i2, d4_ij, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x39



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x39(s_a, i_a, s_i, i_i, s_m, i_m, V2_, Hdiag_, D4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_i, s_i
integer, intent(in) :: i_m, s_m
integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_k, i_k, s_o1, i_o1, s_o2, i_o2
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D4(i,m,k,o1,m,i,o2,k) V2(a,o1,o2,a) 
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),IEOR(s_k,s_o1)) == IEOR(IEOR(s_m,s_i),IEOR(s_o2,s_k)) .and. &
IEOR(s_a,s_o1) == IEOR(s_o2,s_a)) then
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D4_(s_k, s_o2, s_i, s_m, s_o1, s_k)%array(i_k, i_o2, i_i, i_m, i_o1, i_k) & 
  * V2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x39



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x40(sa, ia, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x40(sa, ia, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x40



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x40(s_a, i_a, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_o2, i_o2, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D3(i,i,k,o1,o2,k) V2(a,o1,o2,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_i),s_k) == IEOR(IEOR(s_o1,s_o2),s_k) .and. &
IEOR(s_a,s_o1) == IEOR(s_o2,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D3_(s_k, s_o2, s_o1, s_k, s_i, s_i)%array(i_k, i_o2, i_o1, i_k, i_i, i_i) & 
  * V2_(s_a, s_o2, s_o1)%array(i_a, i_o2, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x40



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x41(sa, ia, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x41(sa, ia, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x41



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x41(s_a, i_a, V2_, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_o1, i_o1, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (    2.00000000) D3(i,k,k,o1,m,i) V2(a,m,o1,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_o1 = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_k),s_k) == IEOR(IEOR(s_o1,s_m),s_i) .and. &
IEOR(s_a,s_m) == IEOR(s_o1,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 2.00000000d+00 & 
  * D3_(s_i, s_m, s_o1, s_k, s_k, s_i)%array(i_i, i_m, i_o1, i_k, i_k, i_i) & 
  * V2_(s_a, s_o1, s_m)%array(i_a, i_o1, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x41



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x42(sa, ia, V2, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x42(sa, ia, h2_i, av2_i2, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x42



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x42(s_a, i_a, V2_, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) D2(i,k,k,i) V2(a,m,m,a) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_k) == IEOR(s_k,s_i) .and. &
IEOR(s_a,s_m) == IEOR(s_m,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * D2_(s_i, s_k, s_k, s_i)%array(i_i, i_k, i_k, i_i) & 
  * V2_(s_a, s_m, s_m)%array(i_a, i_m, i_m)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x42



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x43(sa, ia, Ecas, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Ecas
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x43(sa, ia, Ecas, av2_i2, d3, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x43



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x43(s_a, i_a, Ecas, Hdiag_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Ecas
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_m, i_m, s_k, i_k
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) Ecas D3(i,m,k,k,m,i) 
do s_i = 0, nir-1
do s_m = 0, nir-1
do s_k = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(IEOR(s_i,s_m),s_k) == IEOR(IEOR(s_k,s_m),s_i)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * Ecas & 
  * D3_(s_i, s_m, s_k, s_k, s_m, s_i)%array(i_i, i_m, i_k, i_k, i_m, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x43



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_ooov_no0_x44(sa, ia, Ecas, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: Ecas
real(kind=8), intent(inout) :: Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2_2(sa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_ooov_no0_x44(sa, ia, Ecas, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i2)

end subroutine g_if_diag_ooov_no0_x44



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_diag_ooov_no0_x44(s_a, i_a, Ecas, Hdiag_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Ecas
! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_k, i_k, s_m, i_m
! Hdiag(i,k,m,a)  <-- 
! (    1.00000000) Ecas D2(i,i,k,k) 
do s_i = 0, nir-1
do s_k = 0, nir-1
do s_m = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_m,s_a) .and. & 
IEOR(s_i,s_i) == IEOR(s_k,s_k)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_m = psym(I_BEGIN, I_O, s_m), psym(I_END, I_O, s_m)
Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) =  &
    Hdiag_(s_m, s_k, s_i)%array(i_m, i_k, i_i) &
  + 1.00000000d+00 & 
  * Ecas & 
  * D2_(s_k, s_k, s_i, s_i)%array(i_k, i_k, i_i, i_i)
end do ! Irrep Loop
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_diag_ooov_no0_x44

