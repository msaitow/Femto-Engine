#include "../f_ct.fh"


!    o__ __o__/_                            o                         
!   <|    v                                <|>                        
!   < >                                    < >                        
!    |         o__  __o   \o__ __o__ __o    |        o__ __o         
!    o__/_    /v      |>   |     |     |>   o__/_   /v     v\        
!    |       />      //   / \   / \   / \   |      />       <\    
!   <o>      \o    o/     \o/   \o/   \o/   |      \         /   
!    |        v\  /v __o   |     |     |    o       o       o        
!   / \        <\/> __/>  / \   / \   / \   <\__    <\__ __/>  



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ooov_no0_x0(sv1, iv1, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv1, iv1
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_g_ooov_no0_x0(sv1, iv1, av2_i, h1, xaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h1)
deallocate(xaaaa)

end subroutine g_if_sigma_g_ooov_no0_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_ooov_no0_x0(s_v1, i_v1, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o2, i_o2, s_o1, i_o1, s_o4, i_o4, s_o3, i_o3
! X(o2,o1,o4,o3)  <-- 
! (    1.00000000)  T2(o2,o1,o4,v1) h(o3,v1) 
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o2,s_o1) == IEOR(s_o4,s_o3) .and. & 
IEOR(s_o2,s_o1) == IEOR(s_o4,s_v1) .and. &
IEOR(s_o3,s_v1) == 0) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o4, s_o1, s_o2)%array(i_o3, i_o4, i_o1, i_o2) =  &
    X_(s_o3, s_o4, s_o1, s_o2)%array(i_o3, i_o4, i_o1, i_o2) &
  + 1.00000000d+00 & 
  * T2_(s_o4, s_o1, s_o2)%array(i_o4, i_o1, i_o2) & 
  * h_(s_v1, s_o3)%array(i_v1, i_o3)
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

end subroutine g_sigma_g_ooov_no0_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ooov_no1_x0(X, S0, nir, nsym, psym)

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
call g_sigma_g_ooov_no1_x0(xaaaa, S0, d2, nir, nsym, psym)

deallocate(xaaaa)

end subroutine g_if_sigma_g_ooov_no1_x0



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_ooov_no1_x0(X_, S0_, D2_, nir, nsym, psym)

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
! (    1.00000000) D2(o1,o3,o2,o4) X(o2,o1,o4,o3) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_o2,s_o1) == IEOR(s_o4,s_o3)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
S0_ = S0_ &
  + 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * X_(s_o3, s_o4, s_o1, s_o2)%array(i_o3, i_o4, i_o1, i_o2)
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

end subroutine g_sigma_g_ooov_no1_x0



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ooov_y1(sc1, ic1, V2, Y0, nir, nsym, psym)

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
call set_symblock_Yav(sleft2, Y0, nir, nsym, psym) ! -> Yav (allocate) 
call g_sigma_g_ooov_y1(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_g_ooov_y1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_ooov_y1(s_c1, i_c1, V2_, Y0_, nir, nsym, psym)

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

integer :: s_o2, i_o2, s_v1, i_v1
! Y0(o2,v1)  <-- 
! (    1.00000000)  V2(c1,c1,o2,v1) 
do s_o2 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_o2,s_v1) == 0 .and. & 
IEOR(s_c1,s_c1) == IEOR(s_o2,s_v1)) then
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y0_(s_v1, s_o2)%array(i_v1, i_o2) =  &
    Y0_(s_v1, s_o2)%array(i_v1, i_o2) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_o2, s_c1)%array(i_v1, i_o2, i_c1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_g_ooov_y1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ooov_no0_x1(sv1, iv1, T2, Y0, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y0(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yav(sleft2, Y0, nir, nsym, psym) ! -> Yav (allocate) 
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_g_ooov_no0_x1(sv1, iv1, av2_i, Yav, xaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yav)
deallocate(xaaaa)

end subroutine g_if_sigma_g_ooov_no0_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_ooov_no0_x1(s_v1, i_v1, T2_, Y0_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y0_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o3, i_o3, s_o4, i_o4, s_o1, i_o1, s_o2, i_o2
! X(o3,o4,o1,o2)  <-- 
! (    1.00000000)  T2(o3,o4,o1,v1) Y0(o2,v1) 
do s_o3 = 0, nir-1
do s_o4 = 0, nir-1
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
if( &
IEOR(s_o3,s_o4) == IEOR(s_o1,s_o2) .and. & 
IEOR(s_o3,s_o4) == IEOR(s_o1,s_v1) .and. &
IEOR(s_o2,s_v1) == 0) then
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
X_(s_o2, s_o1, s_o4, s_o3)%array(i_o2, i_o1, i_o4, i_o3) =  &
    X_(s_o2, s_o1, s_o4, s_o3)%array(i_o2, i_o1, i_o4, i_o3) &
  + 1.00000000d+00 & 
  * T2_(s_o1, s_o4, s_o3)%array(i_o1, i_o4, i_o3) & 
  * Y0_(s_v1, s_o2)%array(i_v1, i_o2)
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

end subroutine g_sigma_g_ooov_no0_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ooov_no1_x1(X, S0, nir, nsym, psym)

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
call g_sigma_g_ooov_no1_x1(xaaaa, S0, d2, nir, nsym, psym)

deallocate(xaaaa)

end subroutine g_if_sigma_g_ooov_no1_x1



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_ooov_no1_x1(X_, S0_, D2_, nir, nsym, psym)

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
! (    2.00000000) D2(o1,o3,o2,o4) X(o3,o4,o1,o2) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_o3,s_o4) == IEOR(s_o1,s_o2)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
S0_ = S0_ &
  + 2.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * X_(s_o2, s_o1, s_o4, s_o3)%array(i_o2, i_o1, i_o4, i_o3)
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

end subroutine g_sigma_g_ooov_no1_x1



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ooov_y2(sc1, ic1, V2, Y1, nir, nsym, psym)

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
call set_symblock_Yav(sleft2, Y1, nir, nsym, psym) ! -> Yav (allocate) 
call g_sigma_g_ooov_y2(sc1, ic1, h2_i, Yav, nir, nsym, psym)

deallocate(h2_i)
deallocate(Yav)

end subroutine g_if_sigma_g_ooov_y2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_ooov_y2(s_c1, i_c1, V2_, Y1_, nir, nsym, psym)

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

integer :: s_o1, i_o1, s_v1, i_v1
! Y1(o1,v1)  <-- 
! (    1.00000000)  V2(c1,o1,c1,v1) 
do s_o1 = 0, nir-1
do s_v1 = 0, nir-1
if( &
IEOR(s_o1,s_v1) == 0 .and. & 
IEOR(s_c1,s_o1) == IEOR(s_c1,s_v1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_v1 = psym(I_BEGIN, I_V, s_v1), psym(I_END, I_V, s_v1)
Y1_(s_v1, s_o1)%array(i_v1, i_o1) =  &
    Y1_(s_v1, s_o1)%array(i_v1, i_o1) &
  + 1.00000000d+00 & 
  * V2_(s_v1, s_c1, s_o1)%array(i_v1, i_c1, i_o1)
end do ! Irrep Loop
end do ! Irrep Loop
end if ! Irrep Cond
end do ! Orbital Loop
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_g_ooov_y2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ooov_no0_x2(sv1, iv1, T2, Y1, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv1, iv1
real(kind=8), intent(inout) :: T2(*), Y1(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft2 = 0
call set_symblock_Yav(sleft2, Y1, nir, nsym, psym) ! -> Yav (allocate) 
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_g_ooov_no0_x2(sv1, iv1, av2_i, Yav, xaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(Yav)
deallocate(xaaaa)

end subroutine g_if_sigma_g_ooov_no0_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_ooov_no0_x2(s_v1, i_v1, T2_, Y1_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: Y1_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o4, i_o4, s_o3, i_o3, s_o2, i_o2, s_o1, i_o1
! X(o4,o3,o2,o1)  <-- 
! (    1.00000000)  T2(o4,o3,o2,v1) Y1(o1,v1) 
do s_o4 = 0, nir-1
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o1 = 0, nir-1
if( &
IEOR(s_o4,s_o3) == IEOR(s_o2,s_o1) .and. & 
IEOR(s_o4,s_o3) == IEOR(s_o2,s_v1) .and. &
IEOR(s_o1,s_v1) == 0) then
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
X_(s_o1, s_o2, s_o3, s_o4)%array(i_o1, i_o2, i_o3, i_o4) =  &
    X_(s_o1, s_o2, s_o3, s_o4)%array(i_o1, i_o2, i_o3, i_o4) &
  + 1.00000000d+00 & 
  * T2_(s_o2, s_o3, s_o4)%array(i_o2, i_o3, i_o4) & 
  * Y1_(s_v1, s_o1)%array(i_v1, i_o1)
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

end subroutine g_sigma_g_ooov_no0_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ooov_no1_x2(X, S0, nir, nsym, psym)

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
call g_sigma_g_ooov_no1_x2(xaaaa, S0, d2, nir, nsym, psym)

deallocate(xaaaa)

end subroutine g_if_sigma_g_ooov_no1_x2



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_ooov_no1_x2(X_, S0_, D2_, nir, nsym, psym)

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
! (   -1.00000000) D2(o1,o3,o2,o4) X(o4,o3,o2,o1) 
do s_o1 = 0, nir-1
do s_o3 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
if( &
IEOR(s_o1,s_o3) == IEOR(s_o2,s_o4) .and. &
IEOR(s_o4,s_o3) == IEOR(s_o2,s_o1)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
S0_ = S0_ &
  - 1.00000000d+00 & 
  * D2_(s_o4, s_o2, s_o3, s_o1)%array(i_o4, i_o2, i_o3, i_o1) & 
  * X_(s_o1, s_o2, s_o3, s_o4)%array(i_o1, i_o2, i_o3, i_o4)
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

end subroutine g_sigma_g_ooov_no1_x2



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ooov_no0_x3(sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_xaaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaaa (allocate) 
call g_sigma_g_ooov_no0_x3(sv1, iv1, av2_i, h2_i, xaaaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaaaa)

end subroutine g_if_sigma_g_ooov_no0_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_ooov_no0_x3(s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock6), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_o4, i_o4, s_o5, i_o5, s_o3, i_o3
integer :: s_o6, i_o6
! X(o1,o2,o4,o5,o3,o6)  <-- 
! (    1.00000000)  T2(o1,o2,o4,v1) V2(v1,o5,o3,o6) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_o4 = 0, nir-1
do s_o5 = 0, nir-1
do s_o3 = 0, nir-1
do s_o6 = 0, nir-1
if( &
IEOR(IEOR(s_o1,s_o2),s_o4) == IEOR(IEOR(s_o5,s_o3),s_o6) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_o4,s_v1) .and. &
IEOR(s_v1,s_o5) == IEOR(s_o3,s_o6)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o6 = psym(I_BEGIN, I_O, s_o6), psym(I_END, I_O, s_o6)
X_(s_o6, s_o3, s_o5, s_o4, s_o2, s_o1)%array(i_o6, i_o3, i_o5, i_o4, i_o2, i_o1) =  &
    X_(s_o6, s_o3, s_o5, s_o4, s_o2, s_o1)%array(i_o6, i_o3, i_o5, i_o4, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_o4, s_o2, s_o1)%array(i_o4, i_o2, i_o1) & 
  * V2_(s_o6, s_o3, s_o5)%array(i_o6, i_o3, i_o5)
end do ! Irrep Loop
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
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_g_ooov_no0_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ooov_no1_x3(X, S0, nir, nsym, psym)

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
call set_symblock_xaaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaaa (allocate) 
call g_sigma_g_ooov_no1_x3(xaaaaaa, S0, d3, nir, nsym, psym)

deallocate(xaaaaaa)

end subroutine g_if_sigma_g_ooov_no1_x3



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_ooov_no1_x3(X_, S0_, D3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
real(kind=8)                   :: S0_
type(symblock6), intent(inout) :: D3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o4, i_o4, s_o2, i_o2, s_o5, i_o5, s_o3, i_o3
integer :: s_o6, i_o6
! S0()  <-- 
! (    1.00000000) D3(o1,o4,o2,o5,o3,o6) X(o1,o2,o4,o5,o3,o6) 
do s_o1 = 0, nir-1
do s_o4 = 0, nir-1
do s_o2 = 0, nir-1
do s_o5 = 0, nir-1
do s_o3 = 0, nir-1
do s_o6 = 0, nir-1
if( &
IEOR(IEOR(s_o1,s_o4),s_o2) == IEOR(IEOR(s_o5,s_o3),s_o6) .and. &
IEOR(IEOR(s_o1,s_o2),s_o4) == IEOR(IEOR(s_o5,s_o3),s_o6)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
do i_o6 = psym(I_BEGIN, I_O, s_o6), psym(I_END, I_O, s_o6)
S0_ = S0_ &
  + 1.00000000d+00 & 
  * D3_(s_o6, s_o3, s_o5, s_o2, s_o4, s_o1)%array(i_o6, i_o3, i_o5, i_o2, i_o4, i_o1) & 
  * X_(s_o6, s_o3, s_o5, s_o4, s_o2, s_o1)%array(i_o6, i_o3, i_o5, i_o4, i_o2, i_o1)
end do ! Irrep Loop
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
end do ! Orbital Loop
! FEMTO END  ****************************************************************

end subroutine g_sigma_g_ooov_no1_x3



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ooov_no0_x4(sv1, iv1, T2, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sv1, iv1
real(kind=8), intent(inout) :: T2(*), V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft
integer :: sleft2

call set_symblock_av2(sv1, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h2(sv1, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_g_ooov_no0_x4(sv1, iv1, av2_i, h2_i, xaaaa, nir, nsym, psym)

deallocate(av2_i)
deallocate(h2_i)
deallocate(xaaaa)

end subroutine g_if_sigma_g_ooov_no0_x4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_ooov_no0_x4(s_v1, i_v1, T2_, V2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_v1, s_v1
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_o1, i_o1, s_o2, i_o2, s_o5, i_o5, s_o4, i_o4, s_o3, i_o3
! X(o1,o2,o4,o3)  <-- 
! (    1.00000000)  T2(o1,o2,o5,v1) V2(v1,o4,o3,o5) 
do s_o1 = 0, nir-1
do s_o2 = 0, nir-1
do s_o5 = 0, nir-1
do s_o4 = 0, nir-1
do s_o3 = 0, nir-1
if( &
IEOR(s_o1,s_o2) == IEOR(s_o4,s_o3) .and. & 
IEOR(s_o1,s_o2) == IEOR(s_o5,s_v1) .and. &
IEOR(s_v1,s_o4) == IEOR(s_o3,s_o5)) then
do i_o1 = psym(I_BEGIN, I_O, s_o1), psym(I_END, I_O, s_o1)
do i_o2 = psym(I_BEGIN, I_O, s_o2), psym(I_END, I_O, s_o2)
do i_o5 = psym(I_BEGIN, I_O, s_o5), psym(I_END, I_O, s_o5)
do i_o4 = psym(I_BEGIN, I_O, s_o4), psym(I_END, I_O, s_o4)
do i_o3 = psym(I_BEGIN, I_O, s_o3), psym(I_END, I_O, s_o3)
X_(s_o3, s_o4, s_o2, s_o1)%array(i_o3, i_o4, i_o2, i_o1) =  &
    X_(s_o3, s_o4, s_o2, s_o1)%array(i_o3, i_o4, i_o2, i_o1) &
  + 1.00000000d+00 & 
  * T2_(s_o5, s_o2, s_o1)%array(i_o5, i_o2, i_o1) & 
  * V2_(s_o5, s_o3, s_o4)%array(i_o5, i_o3, i_o4)
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

end subroutine g_sigma_g_ooov_no0_x4



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_g_ooov_no1_x4(X, S0, nir, nsym, psym)

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
call g_sigma_g_ooov_no1_x4(xaaaa, S0, d2, nir, nsym, psym)

deallocate(xaaaa)

end subroutine g_if_sigma_g_ooov_no1_x4



! ---------------------------- Parameters used -------------------------------
!                                                                             
! Whether the LHS is a BareAmpPack ....... No
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                             
!-----------------------------------------------------------------------------

subroutine g_sigma_g_ooov_no1_x4(X_, S0_, D2_, nir, nsym, psym)

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

end subroutine g_sigma_g_ooov_no1_x4

