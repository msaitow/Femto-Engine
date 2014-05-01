

#ifndef C_DIAG_CCVV_H
#define C_DIAG_CCVV_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//  8888888888                     888                  
//  888                            888                  
//  888                            888                  
//  8888888  .d88b.  88888b.d88b.  888888  .d88b.       
//  888     d8P  Y8b 888 "888 "88b 888    d88""88b  
//  888     88888888 888  888  888 888    888  888      
//  888     Y8b.     888  888  888 Y88b.  Y88..88P      
//  888      "Y8888  888  888  888  "Y888  "Y88P"   

void FC_FUNC(g_if_diag_ccvv_y0, G_IF_DIAG_CCVV_Y0)
  (const double * const h, const double * const Y0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x0, G_IF_DIAG_CCVV_NO0_X0)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y0, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y1, G_IF_DIAG_CCVV_Y1)
  (const double * const h, const double * const Y1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x1, G_IF_DIAG_CCVV_NO0_X1)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y1, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x2, G_IF_DIAG_CCVV_NO0_X2)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x3, G_IF_DIAG_CCVV_NO0_X3)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x4, G_IF_DIAG_CCVV_NO0_X4)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y5, G_IF_DIAG_CCVV_Y5)
  (const double * const h, const double * const Y2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x5, G_IF_DIAG_CCVV_NO0_X5)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y2, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y6, G_IF_DIAG_CCVV_Y6)
  (const double * const h, const double * const Y3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x6, G_IF_DIAG_CCVV_NO0_X6)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y3, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x7, G_IF_DIAG_CCVV_NO0_X7)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x8, G_IF_DIAG_CCVV_NO0_X8)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x9, G_IF_DIAG_CCVV_NO0_X9)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x10, G_IF_DIAG_CCVV_NO0_X10)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x11, G_IF_DIAG_CCVV_NO0_X11)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x12, G_IF_DIAG_CCVV_NO0_X12)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x13, G_IF_DIAG_CCVV_NO0_X13)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x14, G_IF_DIAG_CCVV_NO0_X14)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x15, G_IF_DIAG_CCVV_NO0_X15)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x16, G_IF_DIAG_CCVV_NO0_X16)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x17, G_IF_DIAG_CCVV_NO0_X17)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x18, G_IF_DIAG_CCVV_NO0_X18)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y19, G_IF_DIAG_CCVV_Y19)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x19, G_IF_DIAG_CCVV_NO0_X19)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y4, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y20, G_IF_DIAG_CCVV_Y20)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x20, G_IF_DIAG_CCVV_NO0_X20)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y5, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y21, G_IF_DIAG_CCVV_Y21)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x21, G_IF_DIAG_CCVV_NO0_X21)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y6, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y22, G_IF_DIAG_CCVV_Y22)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y7, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x22, G_IF_DIAG_CCVV_NO0_X22)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y7, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y23, G_IF_DIAG_CCVV_Y23)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y8, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x23, G_IF_DIAG_CCVV_NO0_X23)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y8, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y24, G_IF_DIAG_CCVV_Y24)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y9, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x24, G_IF_DIAG_CCVV_NO0_X24)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y9, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x25, G_IF_DIAG_CCVV_NO0_X25)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x26, G_IF_DIAG_CCVV_NO0_X26)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x27, G_IF_DIAG_CCVV_NO0_X27)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const h, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y28, G_IF_DIAG_CCVV_Y28)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y10, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x28, G_IF_DIAG_CCVV_NO0_X28)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y10, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y29, G_IF_DIAG_CCVV_Y29)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y11, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x29, G_IF_DIAG_CCVV_NO0_X29)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y11, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y30, G_IF_DIAG_CCVV_Y30)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y12, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x30, G_IF_DIAG_CCVV_NO0_X30)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y12, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y31, G_IF_DIAG_CCVV_Y31)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y13, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x31, G_IF_DIAG_CCVV_NO0_X31)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y13, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y32, G_IF_DIAG_CCVV_Y32)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y14, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x32, G_IF_DIAG_CCVV_NO0_X32)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y14, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y33, G_IF_DIAG_CCVV_Y33)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y15, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x33, G_IF_DIAG_CCVV_NO0_X33)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y15, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y34, G_IF_DIAG_CCVV_Y34)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y16, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x34, G_IF_DIAG_CCVV_NO0_X34)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y16, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y35, G_IF_DIAG_CCVV_Y35)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y17, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x35, G_IF_DIAG_CCVV_NO0_X35)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y17, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x36, G_IF_DIAG_CCVV_NO0_X36)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x37, G_IF_DIAG_CCVV_NO0_X37)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y38, G_IF_DIAG_CCVV_Y38)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y18, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x38, G_IF_DIAG_CCVV_NO0_X38)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y18, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x39, G_IF_DIAG_CCVV_NO0_X39)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x40, G_IF_DIAG_CCVV_NO0_X40)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y41, G_IF_DIAG_CCVV_Y41)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y19, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x41, G_IF_DIAG_CCVV_NO0_X41)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y19, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x42, G_IF_DIAG_CCVV_NO0_X42)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x43, G_IF_DIAG_CCVV_NO0_X43)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x44, G_IF_DIAG_CCVV_NO0_X44)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y45, G_IF_DIAG_CCVV_Y45)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y20, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x45, G_IF_DIAG_CCVV_NO0_X45)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y20, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x46, G_IF_DIAG_CCVV_NO0_X46)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x47, G_IF_DIAG_CCVV_NO0_X47)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x48, G_IF_DIAG_CCVV_NO0_X48)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y49, G_IF_DIAG_CCVV_Y49)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y21, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x49, G_IF_DIAG_CCVV_NO0_X49)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y21, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x50, G_IF_DIAG_CCVV_NO0_X50)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x51, G_IF_DIAG_CCVV_NO0_X51)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y52, G_IF_DIAG_CCVV_Y52)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y22, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x52, G_IF_DIAG_CCVV_NO0_X52)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y22, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y53, G_IF_DIAG_CCVV_Y53)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y23, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x53, G_IF_DIAG_CCVV_NO0_X53)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y23, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y54, G_IF_DIAG_CCVV_Y54)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y24, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x54, G_IF_DIAG_CCVV_NO0_X54)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y24, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y55, G_IF_DIAG_CCVV_Y55)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y25, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x55, G_IF_DIAG_CCVV_NO0_X55)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y25, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y56, G_IF_DIAG_CCVV_Y56)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y26, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x56, G_IF_DIAG_CCVV_NO0_X56)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y26, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y57, G_IF_DIAG_CCVV_Y57)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y27, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x57, G_IF_DIAG_CCVV_NO0_X57)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y27, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y58, G_IF_DIAG_CCVV_Y58)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y28, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x58, G_IF_DIAG_CCVV_NO0_X58)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y28, 
   const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y59, G_IF_DIAG_CCVV_Y59)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y29, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x59, G_IF_DIAG_CCVV_NO0_X59)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y29, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y60, G_IF_DIAG_CCVV_Y60)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y30, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x60, G_IF_DIAG_CCVV_NO0_X60)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y30, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x61, G_IF_DIAG_CCVV_NO0_X61)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x62, G_IF_DIAG_CCVV_NO0_X62)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y63, G_IF_DIAG_CCVV_Y63)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y31, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x63, G_IF_DIAG_CCVV_NO0_X63)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y31, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y64, G_IF_DIAG_CCVV_Y64)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y32, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x64, G_IF_DIAG_CCVV_NO0_X64)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y32, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y65, G_IF_DIAG_CCVV_Y65)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y33, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x65, G_IF_DIAG_CCVV_NO0_X65)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y33, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x66, G_IF_DIAG_CCVV_NO0_X66)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x67, G_IF_DIAG_CCVV_NO0_X67)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x68, G_IF_DIAG_CCVV_NO0_X68)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sy, const FC_INT &iy, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y69, G_IF_DIAG_CCVV_Y69)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y34, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x69, G_IF_DIAG_CCVV_NO0_X69)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y34, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y70, G_IF_DIAG_CCVV_Y70)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y35, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x70, G_IF_DIAG_CCVV_NO0_X70)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y35, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x71, G_IF_DIAG_CCVV_NO0_X71)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x72, G_IF_DIAG_CCVV_NO0_X72)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sy, const FC_INT &iy, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x73, G_IF_DIAG_CCVV_NO0_X73)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x74, G_IF_DIAG_CCVV_NO0_X74)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x75, G_IF_DIAG_CCVV_NO0_X75)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x76, G_IF_DIAG_CCVV_NO0_X76)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sy, const FC_INT &iy, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y77, G_IF_DIAG_CCVV_Y77)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y36, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x77, G_IF_DIAG_CCVV_NO0_X77)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y36, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x78, G_IF_DIAG_CCVV_NO0_X78)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x79, G_IF_DIAG_CCVV_NO0_X79)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x80, G_IF_DIAG_CCVV_NO0_X80)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x81, G_IF_DIAG_CCVV_NO0_X81)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y82, G_IF_DIAG_CCVV_Y82)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y37, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x82, G_IF_DIAG_CCVV_NO0_X82)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y37, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x83, G_IF_DIAG_CCVV_NO0_X83)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x84, G_IF_DIAG_CCVV_NO0_X84)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x85, G_IF_DIAG_CCVV_NO0_X85)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x86, G_IF_DIAG_CCVV_NO0_X86)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x87, G_IF_DIAG_CCVV_NO0_X87)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sy, const FC_INT &iy, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x88, G_IF_DIAG_CCVV_NO0_X88)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y89, G_IF_DIAG_CCVV_Y89)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y38, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x89, G_IF_DIAG_CCVV_NO0_X89)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y38, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x90, G_IF_DIAG_CCVV_NO0_X90)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y91, G_IF_DIAG_CCVV_Y91)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y39, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x91, G_IF_DIAG_CCVV_NO0_X91)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y39, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y92, G_IF_DIAG_CCVV_Y92)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y40, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x92, G_IF_DIAG_CCVV_NO0_X92)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y40, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y93, G_IF_DIAG_CCVV_Y93)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y41, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x93, G_IF_DIAG_CCVV_NO0_X93)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y41, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x94, G_IF_DIAG_CCVV_NO0_X94)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x95, G_IF_DIAG_CCVV_NO0_X95)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x96, G_IF_DIAG_CCVV_NO0_X96)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x97, G_IF_DIAG_CCVV_NO0_X97)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x98, G_IF_DIAG_CCVV_NO0_X98)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x99, G_IF_DIAG_CCVV_NO0_X99)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x100, G_IF_DIAG_CCVV_NO0_X100)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x101, G_IF_DIAG_CCVV_NO0_X101)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y102, G_IF_DIAG_CCVV_Y102)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y42, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x102, G_IF_DIAG_CCVV_NO0_X102)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y42, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x103, G_IF_DIAG_CCVV_NO0_X103)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y104, G_IF_DIAG_CCVV_Y104)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y43, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x104, G_IF_DIAG_CCVV_NO0_X104)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y43, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x105, G_IF_DIAG_CCVV_NO0_X105)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x106, G_IF_DIAG_CCVV_NO0_X106)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x107, G_IF_DIAG_CCVV_NO0_X107)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sy, const FC_INT &iy, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x108, G_IF_DIAG_CCVV_NO0_X108)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y109, G_IF_DIAG_CCVV_Y109)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y44, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x109, G_IF_DIAG_CCVV_NO0_X109)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y44, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x110, G_IF_DIAG_CCVV_NO0_X110)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y111, G_IF_DIAG_CCVV_Y111)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y45, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x111, G_IF_DIAG_CCVV_NO0_X111)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y45, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y112, G_IF_DIAG_CCVV_Y112)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y46, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x112, G_IF_DIAG_CCVV_NO0_X112)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y46, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y113, G_IF_DIAG_CCVV_Y113)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y47, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x113, G_IF_DIAG_CCVV_NO0_X113)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y47, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x114, G_IF_DIAG_CCVV_NO0_X114)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x115, G_IF_DIAG_CCVV_NO0_X115)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sy, const FC_INT &iy, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x116, G_IF_DIAG_CCVV_NO0_X116)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x117, G_IF_DIAG_CCVV_NO0_X117)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &so1, const FC_INT &io1, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x118, G_IF_DIAG_CCVV_NO0_X118)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &so1, const FC_INT &io1, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x119, G_IF_DIAG_CCVV_NO0_X119)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &so1, const FC_INT &io1, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x120, G_IF_DIAG_CCVV_NO0_X120)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &so1, const FC_INT &io1, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x121, G_IF_DIAG_CCVV_NO0_X121)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x122, G_IF_DIAG_CCVV_NO0_X122)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x123, G_IF_DIAG_CCVV_NO0_X123)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x124, G_IF_DIAG_CCVV_NO0_X124)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x125, G_IF_DIAG_CCVV_NO0_X125)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x126, G_IF_DIAG_CCVV_NO0_X126)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x127, G_IF_DIAG_CCVV_NO0_X127)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y128, G_IF_DIAG_CCVV_Y128)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y48, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x128, G_IF_DIAG_CCVV_NO0_X128)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y48, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x129, G_IF_DIAG_CCVV_NO0_X129)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y130, G_IF_DIAG_CCVV_Y130)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y49, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x130, G_IF_DIAG_CCVV_NO0_X130)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y49, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y131, G_IF_DIAG_CCVV_Y131)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y50, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x131, G_IF_DIAG_CCVV_NO0_X131)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y50, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_y132, G_IF_DIAG_CCVV_Y132)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y51, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x132, G_IF_DIAG_CCVV_NO0_X132)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Y51, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x133, G_IF_DIAG_CCVV_NO0_X133)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x134, G_IF_DIAG_CCVV_NO0_X134)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x135, G_IF_DIAG_CCVV_NO0_X135)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x136, G_IF_DIAG_CCVV_NO0_X136)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x137, G_IF_DIAG_CCVV_NO0_X137)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x138, G_IF_DIAG_CCVV_NO0_X138)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_diag_ccvv_no0_x139, G_IF_DIAG_CCVV_NO0_X139)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const Hdiag, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

      
 }     
       
       
 #endif
       
       
 