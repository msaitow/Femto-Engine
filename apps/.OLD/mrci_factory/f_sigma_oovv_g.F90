#include "../f_ct.fh"


!  8 8888888888   8 8888888888            ,8.       ,8.    8888888 8888888888 ,o888888o.     
!  8 8888         8 8888                 ,888.     ,888.         8 8888    . 8888     `88.   
!  8 8888         8 8888                .`8888.   .`8888.        8 8888   ,8 8888       `8b  
!  8 8888         8 8888               ,8.`8888. ,8.`8888.       8 8888   88 8888        `8b 
!  8 888888888888 8 888888888888      ,8'8.`8888,8^8.`8888.      8 8888   88 8888         88 
!  8 8888         8 8888             ,8' `8.`8888' `8.`8888.     8 8888   88 8888         88 
!  8 8888         8 8888            ,8'   `8.`88'   `8.`8888.    8 8888   88 8888        ,8P 
!  8 8888         8 8888           ,8'     `8.`'     `8.`8888.   8 8888   `8 8888       ,8P  
!  8 8888         8 8888          ,8'       `8        `8.`8888.  8 8888    ` 8888     ,88'   
!  8 8888         8 888888888888 ,8'         `         `8.`8888. 8 8888       `8888888P'     



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_g_no0_x0(sc, ic, T0, V2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic
real(kind=8), intent(inout) :: T0
real(kind=8), intent(inout) :: V2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sc, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_g_no0_x0(sc, ic, T0, h2_i, av2_i2, d2, nir, nsym, psym)

deallocate(h2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_g_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... Yes
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_oovv_g_no0_x0(s_c, i_c, T0, V2_, S2_, D2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: T0
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: D2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_i, i_i, s_o1, i_o1, s_k, i_k, s_o2, i_o2, s_a, i_a
! S2(i,k,a,c)  <-- 
! (    1.00000000) T0 D2(i,o1,k,o2) V2(c,o2,o1,a) 
do s_i = 0, nir-1
do s_o1 = 0, nir-1
do s_k = 0, nir-1
do s_o2 = 0, nir-1
do s_a = 0, nir-1
if( &
IEOR(s_i,s_k) == IEOR(s_a,s_c) .and. & 
IEOR(s_i,s_o1) == IEOR(s_k,s_o2) .and. &
IEOR(s_c,s_o2) == IEOR(s_o1,s_a)) then
do i_i = psym(I_BEGIN, I_O, s_i), psym(I_END, I_O, s_i)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_k = psym(I_BEGIN, I_O, s_k), psym(I_END, I_O, s_k)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_a = psym(I_BEGIN, I_V, s_a), psym(I_END, I_V, s_a)
S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) =  &
    S2_(s_a, s_k, s_i)%array(i_a, i_k, i_i) &
  + 1.00000000d+00 & 
  * T0 & 
  * D2_(s_o2, s_k, s_o1, s_i)%array(i_o2, i_k, i_o1, i_i) & 
  * V2_(s_a, s_o1, s_o2)%array(i_a, i_o1, i_o2)
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

end subroutine g_sigma_oovv_g_no0_x0

