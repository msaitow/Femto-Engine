

#ifndef C_SIGMA_OOVV_OOVV_H
#define C_SIGMA_OOVV_OOVV_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//       #                #########      #     #   # 
//  ########## ##########         #   #######  #   # 
//      #    #         #          #    # #     #   # 
//      #    #        #   ########     # #     #   # 
//     #     #     # #           #  ##########    #  
//    #   # #       #            #       #       #   
//   #     #         #    ########       #     ##    

void FC_FUNC(g_if_sigma_oovv_oovv_y0, G_IF_SIGMA_OOVV_OOVV_Y0)
  (const double * const h, const double * const Y0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x0, G_IF_SIGMA_OOVV_OOVV_NO0_X0)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y1, G_IF_SIGMA_OOVV_OOVV_Y1)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x1, G_IF_SIGMA_OOVV_OOVV_NO0_X1)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y1, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y2, G_IF_SIGMA_OOVV_OOVV_Y2)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x2, G_IF_SIGMA_OOVV_OOVV_NO0_X2)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const Y2, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y3, G_IF_SIGMA_OOVV_OOVV_Y3)
  (const double * const h, const double * const Y3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x3, G_IF_SIGMA_OOVV_OOVV_NO0_X3)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const Y3, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y4, G_IF_SIGMA_OOVV_OOVV_Y4)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x4, G_IF_SIGMA_OOVV_OOVV_NO0_X4)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const Y4, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y5, G_IF_SIGMA_OOVV_OOVV_Y5)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x5, G_IF_SIGMA_OOVV_OOVV_NO0_X5)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const Y5, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x6, G_IF_SIGMA_OOVV_OOVV_NO0_X6)
  (const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x6, G_IF_SIGMA_OOVV_OOVV_NO1_X6)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x7, G_IF_SIGMA_OOVV_OOVV_NO0_X7)
  (const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x7, G_IF_SIGMA_OOVV_OOVV_NO1_X7)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x8, G_IF_SIGMA_OOVV_OOVV_NO0_X8)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x8, G_IF_SIGMA_OOVV_OOVV_NO1_X8)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x9, G_IF_SIGMA_OOVV_OOVV_NO0_X9)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x9, G_IF_SIGMA_OOVV_OOVV_NO1_X9)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x10, G_IF_SIGMA_OOVV_OOVV_NO0_X10)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x10, G_IF_SIGMA_OOVV_OOVV_NO1_X10)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x11, G_IF_SIGMA_OOVV_OOVV_NO0_X11)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x11, G_IF_SIGMA_OOVV_OOVV_NO1_X11)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y12, G_IF_SIGMA_OOVV_OOVV_Y12)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x12, G_IF_SIGMA_OOVV_OOVV_NO0_X12)
  (const double * const Y6, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x12, G_IF_SIGMA_OOVV_OOVV_NO1_X12)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y13, G_IF_SIGMA_OOVV_OOVV_Y13)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y7, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x13, G_IF_SIGMA_OOVV_OOVV_NO0_X13)
  (const double * const Y7, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x13, G_IF_SIGMA_OOVV_OOVV_NO1_X13)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y14, G_IF_SIGMA_OOVV_OOVV_Y14)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y8, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x14, G_IF_SIGMA_OOVV_OOVV_NO0_X14)
  (const double * const Y8, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x14, G_IF_SIGMA_OOVV_OOVV_NO1_X14)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y15, G_IF_SIGMA_OOVV_OOVV_Y15)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y9, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x15, G_IF_SIGMA_OOVV_OOVV_NO0_X15)
  (const double * const Y9, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x15, G_IF_SIGMA_OOVV_OOVV_NO1_X15)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y16, G_IF_SIGMA_OOVV_OOVV_Y16)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y10, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x16, G_IF_SIGMA_OOVV_OOVV_NO0_X16)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const Y10, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x16, G_IF_SIGMA_OOVV_OOVV_NO1_X16)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y17, G_IF_SIGMA_OOVV_OOVV_Y17)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y11, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x17, G_IF_SIGMA_OOVV_OOVV_NO0_X17)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const Y11, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x17, G_IF_SIGMA_OOVV_OOVV_NO1_X17)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y18, G_IF_SIGMA_OOVV_OOVV_Y18)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y12, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x18, G_IF_SIGMA_OOVV_OOVV_NO0_X18)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const Y12, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x18, G_IF_SIGMA_OOVV_OOVV_NO1_X18)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y19, G_IF_SIGMA_OOVV_OOVV_Y19)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y13, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x19, G_IF_SIGMA_OOVV_OOVV_NO0_X19)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const Y13, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x19, G_IF_SIGMA_OOVV_OOVV_NO1_X19)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y20, G_IF_SIGMA_OOVV_OOVV_Y20)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y14, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x20, G_IF_SIGMA_OOVV_OOVV_NO0_X20)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const Y14, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x20, G_IF_SIGMA_OOVV_OOVV_NO1_X20)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y21, G_IF_SIGMA_OOVV_OOVV_Y21)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y15, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x21, G_IF_SIGMA_OOVV_OOVV_NO0_X21)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const Y15, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x21, G_IF_SIGMA_OOVV_OOVV_NO1_X21)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y22, G_IF_SIGMA_OOVV_OOVV_Y22)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y16, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x22, G_IF_SIGMA_OOVV_OOVV_NO0_X22)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const Y16, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x22, G_IF_SIGMA_OOVV_OOVV_NO1_X22)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_y23, G_IF_SIGMA_OOVV_OOVV_Y23)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y17, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x23, G_IF_SIGMA_OOVV_OOVV_NO0_X23)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const Y17, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x23, G_IF_SIGMA_OOVV_OOVV_NO1_X23)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x24, G_IF_SIGMA_OOVV_OOVV_NO0_X24)
  (const FC_INT &so1, const FC_INT &io1, const FC_INT &so5, const FC_INT &io5, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x24, G_IF_SIGMA_OOVV_OOVV_NO1_X24)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x25, G_IF_SIGMA_OOVV_OOVV_NO0_X25)
  (const FC_INT &so1, const FC_INT &io1, const FC_INT &so5, const FC_INT &io5, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x25, G_IF_SIGMA_OOVV_OOVV_NO1_X25)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x26, G_IF_SIGMA_OOVV_OOVV_NO0_X26)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x26, G_IF_SIGMA_OOVV_OOVV_NO1_X26)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x27, G_IF_SIGMA_OOVV_OOVV_NO0_X27)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x27, G_IF_SIGMA_OOVV_OOVV_NO1_X27)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x28, G_IF_SIGMA_OOVV_OOVV_NO0_X28)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x28, G_IF_SIGMA_OOVV_OOVV_NO1_X28)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x29, G_IF_SIGMA_OOVV_OOVV_NO0_X29)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x29, G_IF_SIGMA_OOVV_OOVV_NO1_X29)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x30, G_IF_SIGMA_OOVV_OOVV_NO0_X30)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x30, G_IF_SIGMA_OOVV_OOVV_NO1_X30)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x31, G_IF_SIGMA_OOVV_OOVV_NO0_X31)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x31, G_IF_SIGMA_OOVV_OOVV_NO1_X31)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x32, G_IF_SIGMA_OOVV_OOVV_NO0_X32)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x32, G_IF_SIGMA_OOVV_OOVV_NO1_X32)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x33, G_IF_SIGMA_OOVV_OOVV_NO0_X33)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x33, G_IF_SIGMA_OOVV_OOVV_NO1_X33)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x34, G_IF_SIGMA_OOVV_OOVV_NO0_X34)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x34, G_IF_SIGMA_OOVV_OOVV_NO1_X34)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x35, G_IF_SIGMA_OOVV_OOVV_NO0_X35)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sv2, const FC_INT &iv2, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x35, G_IF_SIGMA_OOVV_OOVV_NO1_X35)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

      
 }     
       
       
 #endif
       
       
 