

#ifndef C_SIGMA_OOVV_OOOV_H
#define C_SIGMA_OOVV_OOOV_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//  __/\\\\\\\\\\\\\\\____________________________________________________________________                                   
//   _\/\\\///////////_____________________________________________________________________                                             
//    _\/\\\_______________________________________________________/\\\_____________________                                         
//     _\/\\\\\\\\\\\__________/\\\\\\\\______/\\\\\__/\\\\\_____/\\\\\\\\\\\______/\\\\\____ 
//      _\/\\\///////_________/\\\/////\\\___/\\\///\\\\\///\\\__\////\\\////_____/\\\///\\\__               
//       _\/\\\_______________/\\\\\\\\\\\___\/\\\_\//\\\__\/\\\_____\/\\\________/\\\__\//\\\_       
//        _\/\\\______________\//\\///////____\/\\\__\/\\\__\/\\\_____\/\\\_/\\___\//\\\__/\\\__            
//         _\/\\\_______________\//\\\\\\\\\\__\/\\\__\/\\\__\/\\\_____\//\\\\\_____\///\\\\\/___    
//          _\///_________________\//////////___\///___\///___\///_______\/////________\/////_____                                   

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x0, G_IF_SIGMA_OOVV_OOOV_NO0_X0)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x0, G_IF_SIGMA_OOVV_OOOV_NO1_X0)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const h, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x1, G_IF_SIGMA_OOVV_OOOV_NO0_X1)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x1, G_IF_SIGMA_OOVV_OOOV_NO1_X1)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const h, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x2, G_IF_SIGMA_OOVV_OOOV_NO0_X2)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x2, G_IF_SIGMA_OOVV_OOOV_NO1_X2)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const h, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x3, G_IF_SIGMA_OOVV_OOOV_NO0_X3)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x3, G_IF_SIGMA_OOVV_OOOV_NO1_X3)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const h, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_y4, G_IF_SIGMA_OOVV_OOOV_Y4)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x4, G_IF_SIGMA_OOVV_OOOV_NO0_X4)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x4, G_IF_SIGMA_OOVV_OOOV_NO1_X4)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const Y0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_y5, G_IF_SIGMA_OOVV_OOOV_Y5)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x5, G_IF_SIGMA_OOVV_OOOV_NO0_X5)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x5, G_IF_SIGMA_OOVV_OOOV_NO1_X5)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const Y1, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_y6, G_IF_SIGMA_OOVV_OOOV_Y6)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x6, G_IF_SIGMA_OOVV_OOOV_NO0_X6)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x6, G_IF_SIGMA_OOVV_OOOV_NO1_X6)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const Y2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_y7, G_IF_SIGMA_OOVV_OOOV_Y7)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x7, G_IF_SIGMA_OOVV_OOOV_NO0_X7)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x7, G_IF_SIGMA_OOVV_OOOV_NO1_X7)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const Y3, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_y8, G_IF_SIGMA_OOVV_OOOV_Y8)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x8, G_IF_SIGMA_OOVV_OOOV_NO0_X8)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x8, G_IF_SIGMA_OOVV_OOOV_NO1_X8)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const Y4, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_y9, G_IF_SIGMA_OOVV_OOOV_Y9)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x9, G_IF_SIGMA_OOVV_OOOV_NO0_X9)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x9, G_IF_SIGMA_OOVV_OOOV_NO1_X9)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const Y5, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_y10, G_IF_SIGMA_OOVV_OOOV_Y10)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x10, G_IF_SIGMA_OOVV_OOOV_NO0_X10)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x10, G_IF_SIGMA_OOVV_OOOV_NO1_X10)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const Y6, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_y11, G_IF_SIGMA_OOVV_OOOV_Y11)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y7, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x11, G_IF_SIGMA_OOVV_OOOV_NO0_X11)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x11, G_IF_SIGMA_OOVV_OOOV_NO1_X11)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const Y7, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x12, G_IF_SIGMA_OOVV_OOOV_NO0_X12)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &si, const FC_INT &ii, const FC_INT &so4, const FC_INT &io4, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x12, G_IF_SIGMA_OOVV_OOOV_NO1_X12)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, const FC_INT &si, const FC_INT &ii, const FC_INT &so4, const FC_INT &io4, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x13, G_IF_SIGMA_OOVV_OOOV_NO0_X13)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &so2, const FC_INT &io2, const FC_INT &so6, const FC_INT &io6, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x13, G_IF_SIGMA_OOVV_OOOV_NO1_X13)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &so2, const FC_INT &io2, const FC_INT &so6, const FC_INT &io6, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x14, G_IF_SIGMA_OOVV_OOOV_NO0_X14)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x14, G_IF_SIGMA_OOVV_OOOV_NO1_X14)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x15, G_IF_SIGMA_OOVV_OOOV_NO0_X15)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x15, G_IF_SIGMA_OOVV_OOOV_NO1_X15)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x16, G_IF_SIGMA_OOVV_OOOV_NO0_X16)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x16, G_IF_SIGMA_OOVV_OOOV_NO1_X16)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x17, G_IF_SIGMA_OOVV_OOOV_NO0_X17)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x17, G_IF_SIGMA_OOVV_OOOV_NO1_X17)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x18, G_IF_SIGMA_OOVV_OOOV_NO0_X18)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x18, G_IF_SIGMA_OOVV_OOOV_NO1_X18)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x19, G_IF_SIGMA_OOVV_OOOV_NO0_X19)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x19, G_IF_SIGMA_OOVV_OOOV_NO1_X19)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x20, G_IF_SIGMA_OOVV_OOOV_NO0_X20)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x20, G_IF_SIGMA_OOVV_OOOV_NO1_X20)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no0_x21, G_IF_SIGMA_OOVV_OOOV_NO0_X21)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ooov_no1_x21, G_IF_SIGMA_OOVV_OOOV_NO1_X21)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

      
 }     
       
       
 #endif
       
       
 