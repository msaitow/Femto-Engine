

#ifndef C_SIGMA_OOOV_OOOV_H
#define C_SIGMA_OOOV_OOOV_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//      ______                  __           
//     / ____/___   ____ ___   / /_ ____     
//    / /_   / _ \ / __ `__ \ / __// __ \ 
//   / __/  /  __// / / / / // /_ / /_/ /    
//  /_/     \___//_/ /_/ /_/ \__/ \____/  

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x0, G_IF_SIGMA_OOOV_OOOV_NO0_X0)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x0, G_IF_SIGMA_OOOV_OOOV_NO1_X0)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x1, G_IF_SIGMA_OOOV_OOOV_NO0_X1)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x1, G_IF_SIGMA_OOOV_OOOV_NO1_X1)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x2, G_IF_SIGMA_OOOV_OOOV_NO0_X2)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x2, G_IF_SIGMA_OOOV_OOOV_NO1_X2)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x3, G_IF_SIGMA_OOOV_OOOV_NO0_X3)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x3, G_IF_SIGMA_OOOV_OOOV_NO1_X3)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x4, G_IF_SIGMA_OOOV_OOOV_NO0_X4)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x4, G_IF_SIGMA_OOOV_OOOV_NO1_X4)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x5, G_IF_SIGMA_OOOV_OOOV_NO0_X5)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x5, G_IF_SIGMA_OOOV_OOOV_NO1_X5)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x6, G_IF_SIGMA_OOOV_OOOV_NO0_X6)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x6, G_IF_SIGMA_OOOV_OOOV_NO1_X6)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x7, G_IF_SIGMA_OOOV_OOOV_NO0_X7)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x7, G_IF_SIGMA_OOOV_OOOV_NO1_X7)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y8, G_IF_SIGMA_OOOV_OOOV_Y8)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x8, G_IF_SIGMA_OOOV_OOOV_NO0_X8)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const Y0, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x8, G_IF_SIGMA_OOOV_OOOV_NO1_X8)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y9, G_IF_SIGMA_OOOV_OOOV_Y9)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x9, G_IF_SIGMA_OOOV_OOOV_NO0_X9)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const Y1, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x9, G_IF_SIGMA_OOOV_OOOV_NO1_X9)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y10, G_IF_SIGMA_OOOV_OOOV_Y10)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x10, G_IF_SIGMA_OOOV_OOOV_NO0_X10)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const Y2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x10, G_IF_SIGMA_OOOV_OOOV_NO1_X10)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y11, G_IF_SIGMA_OOOV_OOOV_Y11)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x11, G_IF_SIGMA_OOOV_OOOV_NO0_X11)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const Y3, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x11, G_IF_SIGMA_OOOV_OOOV_NO1_X11)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y12, G_IF_SIGMA_OOOV_OOOV_Y12)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x12, G_IF_SIGMA_OOOV_OOOV_NO0_X12)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const Y4, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x12, G_IF_SIGMA_OOOV_OOOV_NO1_X12)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y13, G_IF_SIGMA_OOOV_OOOV_Y13)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x13, G_IF_SIGMA_OOOV_OOOV_NO0_X13)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const Y5, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x13, G_IF_SIGMA_OOOV_OOOV_NO1_X13)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y14, G_IF_SIGMA_OOOV_OOOV_Y14)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x14, G_IF_SIGMA_OOOV_OOOV_NO0_X14)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const Y6, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x14, G_IF_SIGMA_OOOV_OOOV_NO1_X14)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y15, G_IF_SIGMA_OOOV_OOOV_Y15)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y7, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x15, G_IF_SIGMA_OOOV_OOOV_NO0_X15)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const Y7, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x15, G_IF_SIGMA_OOOV_OOOV_NO1_X15)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y16, G_IF_SIGMA_OOOV_OOOV_Y16)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y8, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x16, G_IF_SIGMA_OOOV_OOOV_NO0_X16)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const Y8, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x16, G_IF_SIGMA_OOOV_OOOV_NO1_X16)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y17, G_IF_SIGMA_OOOV_OOOV_Y17)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y9, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x17, G_IF_SIGMA_OOOV_OOOV_NO0_X17)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const Y9, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x17, G_IF_SIGMA_OOOV_OOOV_NO1_X17)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y18, G_IF_SIGMA_OOOV_OOOV_Y18)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y10, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x18, G_IF_SIGMA_OOOV_OOOV_NO0_X18)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const Y10, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x18, G_IF_SIGMA_OOOV_OOOV_NO1_X18)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y19, G_IF_SIGMA_OOOV_OOOV_Y19)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y11, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x19, G_IF_SIGMA_OOOV_OOOV_NO0_X19)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const Y11, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x19, G_IF_SIGMA_OOOV_OOOV_NO1_X19)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y20, G_IF_SIGMA_OOOV_OOOV_Y20)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y12, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x20, G_IF_SIGMA_OOOV_OOOV_NO0_X20)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const Y12, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x20, G_IF_SIGMA_OOOV_OOOV_NO1_X20)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y21, G_IF_SIGMA_OOOV_OOOV_Y21)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y13, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x21, G_IF_SIGMA_OOOV_OOOV_NO0_X21)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const Y13, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x21, G_IF_SIGMA_OOOV_OOOV_NO1_X21)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y22, G_IF_SIGMA_OOOV_OOOV_Y22)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y14, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x22, G_IF_SIGMA_OOOV_OOOV_NO0_X22)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const Y14, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x22, G_IF_SIGMA_OOOV_OOOV_NO1_X22)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_y23, G_IF_SIGMA_OOOV_OOOV_Y23)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y15, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x23, G_IF_SIGMA_OOOV_OOOV_NO0_X23)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const Y15, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x23, G_IF_SIGMA_OOOV_OOOV_NO1_X23)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x24, G_IF_SIGMA_OOOV_OOOV_NO0_X24)
  (const FC_INT &so3, const FC_INT &io3, const FC_INT &so6, const FC_INT &io6, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x24, G_IF_SIGMA_OOOV_OOOV_NO1_X24)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &so3, const FC_INT &io3, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x25, G_IF_SIGMA_OOOV_OOOV_NO0_X25)
  (const FC_INT &sm, const FC_INT &im, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x25, G_IF_SIGMA_OOOV_OOOV_NO1_X25)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sm, const FC_INT &im, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x26, G_IF_SIGMA_OOOV_OOOV_NO0_X26)
  (const FC_INT &sm, const FC_INT &im, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x26, G_IF_SIGMA_OOOV_OOOV_NO1_X26)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sm, const FC_INT &im, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x27, G_IF_SIGMA_OOOV_OOOV_NO0_X27)
  (const FC_INT &so3, const FC_INT &io3, const FC_INT &so6, const FC_INT &io6, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x27, G_IF_SIGMA_OOOV_OOOV_NO1_X27)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &so6, const FC_INT &io6, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x28, G_IF_SIGMA_OOOV_OOOV_NO0_X28)
  (const FC_INT &so1, const FC_INT &io1, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x28, G_IF_SIGMA_OOOV_OOOV_NO1_X28)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &so1, const FC_INT &io1, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x29, G_IF_SIGMA_OOOV_OOOV_NO0_X29)
  (const FC_INT &sm, const FC_INT &im, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x29, G_IF_SIGMA_OOOV_OOOV_NO1_X29)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sm, const FC_INT &im, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x30, G_IF_SIGMA_OOOV_OOOV_NO0_X30)
  (const FC_INT &sk, const FC_INT &ik, const FC_INT &so4, const FC_INT &io4, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x30, G_IF_SIGMA_OOOV_OOOV_NO1_X30)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sk, const FC_INT &ik, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x31, G_IF_SIGMA_OOOV_OOOV_NO0_X31)
  (const FC_INT &so1, const FC_INT &io1, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x31, G_IF_SIGMA_OOOV_OOOV_NO1_X31)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &so1, const FC_INT &io1, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x32, G_IF_SIGMA_OOOV_OOOV_NO0_X32)
  (const FC_INT &sm, const FC_INT &im, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x32, G_IF_SIGMA_OOOV_OOOV_NO1_X32)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sm, const FC_INT &im, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x33, G_IF_SIGMA_OOOV_OOOV_NO0_X33)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &so1, const FC_INT &io1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x33, G_IF_SIGMA_OOOV_OOOV_NO1_X33)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x34, G_IF_SIGMA_OOOV_OOOV_NO0_X34)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &so1, const FC_INT &io1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x34, G_IF_SIGMA_OOOV_OOOV_NO1_X34)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x35, G_IF_SIGMA_OOOV_OOOV_NO0_X35)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x35, G_IF_SIGMA_OOOV_OOOV_NO1_X35)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, const FC_INT &sm, const FC_INT &im, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x36, G_IF_SIGMA_OOOV_OOOV_NO0_X36)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x36, G_IF_SIGMA_OOOV_OOOV_NO1_X36)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x37, G_IF_SIGMA_OOOV_OOOV_NO0_X37)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x37, G_IF_SIGMA_OOOV_OOOV_NO1_X37)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x38, G_IF_SIGMA_OOOV_OOOV_NO0_X38)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x38, G_IF_SIGMA_OOOV_OOOV_NO1_X38)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x39, G_IF_SIGMA_OOOV_OOOV_NO0_X39)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x39, G_IF_SIGMA_OOOV_OOOV_NO1_X39)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x40, G_IF_SIGMA_OOOV_OOOV_NO0_X40)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x40, G_IF_SIGMA_OOOV_OOOV_NO1_X40)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, const FC_INT &sm, const FC_INT &im, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x41, G_IF_SIGMA_OOOV_OOOV_NO0_X41)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x41, G_IF_SIGMA_OOOV_OOOV_NO1_X41)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x42, G_IF_SIGMA_OOOV_OOOV_NO0_X42)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x42, G_IF_SIGMA_OOOV_OOOV_NO1_X42)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x43, G_IF_SIGMA_OOOV_OOOV_NO0_X43)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x43, G_IF_SIGMA_OOOV_OOOV_NO1_X43)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x44, G_IF_SIGMA_OOOV_OOOV_NO0_X44)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x44, G_IF_SIGMA_OOOV_OOOV_NO1_X44)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x45, G_IF_SIGMA_OOOV_OOOV_NO0_X45)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Ecas, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x46, G_IF_SIGMA_OOOV_OOOV_NO0_X46)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Ecas, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

      
 }     
       
       
 #endif
       
       
 