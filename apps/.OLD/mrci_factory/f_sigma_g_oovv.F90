#include "../f_ct.fh"


!                                                              
!   _______________                                  ______    
!  |          |                 .'. .`. `````|`````.~      ~.  
!  |______    |______         .'   `   `.    |    |          | 
!  |          |             .'           `.  |    |          | 
!  |          |___________.'               `.|     `.______.'  
!                                                              



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_oovv_no0_x0(sv2, iv2, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv2, iv2
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv2, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv2, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_g_oovv_no0_x0(sv2, iv2, av2_i, h2_i, xaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaa)

end subroutine g_if_sigma_g_oovv_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_oovv_no0_x0(s_v2, i_v2, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v2, s_v2
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_v1, i_v1, s_o4, i_o4, s_o3, i_o3
! X(o1,o2,o4,o3)  <-- 
! (    1.00000000)  T2(o1,o2,v1,v2) V2(v2,o4,o3,v1) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
do s_o4 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_o4,s_o3) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_v1,s_v2) .and. &
IEOR(s_v2,s_o4) == IEOR(s_o3,s_v1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o4, s_o2, s_o1)%array(i_o3, i_o4, i_o2, i_o1) =  &
    X_(s_o3, s_o4, s_o2, s_o1)%array(i_o3, i_o4, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_v1, s_o2, s_o1)%array(i_v1, i_o2, i_o1) & 
  * V2_(s_v1, s_o3, s_o4)%array(i_v1, i_o3, i_o4)
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

end subroutine g_sigma_g_oovv_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_oovv_no1_x0(X, S0, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: X(*), S0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_g_oovv_no1_x0(xaaaa, S0, d2, nir, nsym, psym)

deallocate(xaaaa)

end subroutine g_if_sigma_g_oovv_no1_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_oovv_no1_x0(X_, S0_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: S0_
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o3, i_o3, s_o2, i_o2, s_o4, i_o4
! S0()  <-- 
! (    1.00000000) D2(o1,o3,o2,o4) X(o1,o2,o4,o3) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_o1,s_o2) == IEOR(s_o4,s_o3)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
S0_ = S0_ &
  + 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * X_(s_o3, s_o4, s_o2, s_o1)%array(i_o3, i_o4, i_o2, i_o1)
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

end subroutine g_sigma_g_oovv_no1_x0

