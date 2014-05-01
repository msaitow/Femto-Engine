

#ifndef C_DIAG_OOVV_H
#define C_DIAG_OOVV_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//  8 8888888888   8 8888888888            ,8.       ,8.    8888888 8888888888 ,o888888o.     
//  8 8888         8 8888                 ,888.     ,888.         8 8888    . 8888     `88.   
//  8 8888         8 8888                .`8888.   .`8888.        8 8888   ,8 8888       `8b  
//  8 8888         8 8888               ,8.`8888. ,8.`8888.       8 8888   88 8888        `8b 
//  8 888888888888 8 888888888888      ,8'8.`8888,8^8.`8888.      8 8888   88 8888         88 
//  8 8888         8 8888             ,8' `8.`8888' `8.`8888.     8 8888   88 8888         88 
//  8 8888         8 8888            ,8'   `8.`88'   `8.`8888.    8 8888   88 8888        ,8P 
//  8 8888         8 8888           ,8'     `8.`'     `8.`8888.   8 8888   `8 8888       ,8P  
//  8 8888         8 8888          ,8'       `8        `8.`8888.  8 8888    ` 8888     ,88'   
//  8 8888         8 888888888888 ,8'         `         `8.`8888. 8 8888       `8888888P'     

void FC_FUNC(g_if_diag_oovv_y0, G_IF_DIAG_OOVV_Y0)
  (const double * const h, const double * const Y0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x0, G_IF_DIAG_OOVV_NO0_X0)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y0, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x1, G_IF_DIAG_OOVV_NO0_X1)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x2, G_IF_DIAG_OOVV_NO0_X2)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y3, G_IF_DIAG_OOVV_Y3)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x3, G_IF_DIAG_OOVV_NO0_X3)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y1, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y4, G_IF_DIAG_OOVV_Y4)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x4, G_IF_DIAG_OOVV_NO0_X4)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y2, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y5, G_IF_DIAG_OOVV_Y5)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x5, G_IF_DIAG_OOVV_NO0_X5)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y3, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y6, G_IF_DIAG_OOVV_Y6)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x6, G_IF_DIAG_OOVV_NO0_X6)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y4, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y7, G_IF_DIAG_OOVV_Y7)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x7, G_IF_DIAG_OOVV_NO0_X7)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y5, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y8, G_IF_DIAG_OOVV_Y8)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x8, G_IF_DIAG_OOVV_NO0_X8)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y6, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x9, G_IF_DIAG_OOVV_NO0_X9)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y10, G_IF_DIAG_OOVV_Y10)
  (const double * const h, const double * const Y7, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x10, G_IF_DIAG_OOVV_NO0_X10)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y7, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x11, G_IF_DIAG_OOVV_NO0_X11)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x12, G_IF_DIAG_OOVV_NO0_X12)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x13, G_IF_DIAG_OOVV_NO0_X13)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y14, G_IF_DIAG_OOVV_Y14)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y8, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x14, G_IF_DIAG_OOVV_NO0_X14)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y8, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y15, G_IF_DIAG_OOVV_Y15)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y9, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x15, G_IF_DIAG_OOVV_NO0_X15)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y9, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y16, G_IF_DIAG_OOVV_Y16)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y10, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x16, G_IF_DIAG_OOVV_NO0_X16)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y10, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y17, G_IF_DIAG_OOVV_Y17)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y11, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x17, G_IF_DIAG_OOVV_NO0_X17)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y11, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x18, G_IF_DIAG_OOVV_NO0_X18)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x19, G_IF_DIAG_OOVV_NO0_X19)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y20, G_IF_DIAG_OOVV_Y20)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y12, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x20, G_IF_DIAG_OOVV_NO0_X20)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y12, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y21, G_IF_DIAG_OOVV_Y21)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y13, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x21, G_IF_DIAG_OOVV_NO0_X21)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y13, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y22, G_IF_DIAG_OOVV_Y22)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y14, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x22, G_IF_DIAG_OOVV_NO0_X22)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y14, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y23, G_IF_DIAG_OOVV_Y23)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y15, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x23, G_IF_DIAG_OOVV_NO0_X23)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y15, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y24, G_IF_DIAG_OOVV_Y24)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y16, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x24, G_IF_DIAG_OOVV_NO0_X24)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y16, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y25, G_IF_DIAG_OOVV_Y25)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y17, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x25, G_IF_DIAG_OOVV_NO0_X25)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y17, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x26, G_IF_DIAG_OOVV_NO0_X26)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &so1, const FC_INT &io1, const FC_INT &so3, const FC_INT &io3, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x27, G_IF_DIAG_OOVV_NO0_X27)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &so1, const FC_INT &io1, const FC_INT &so3, const FC_INT &io3, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x28, G_IF_DIAG_OOVV_NO0_X28)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x29, G_IF_DIAG_OOVV_NO0_X29)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x30, G_IF_DIAG_OOVV_NO0_X30)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x31, G_IF_DIAG_OOVV_NO0_X31)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y32, G_IF_DIAG_OOVV_Y32)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y18, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x32, G_IF_DIAG_OOVV_NO0_X32)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y18, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_y33, G_IF_DIAG_OOVV_Y33)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y19, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x33, G_IF_DIAG_OOVV_NO0_X33)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y19, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x34, G_IF_DIAG_OOVV_NO0_X34)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x35, G_IF_DIAG_OOVV_NO0_X35)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_oovv_no0_x36, G_IF_DIAG_OOVV_NO0_X36)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

      
 }     
       
       
 #endif
       
       
 