

#ifndef C_SIGMA_CCVV_CCOV_H
#define C_SIGMA_CCVV_CCOV_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//  `7MM"""YMM                         mm               
//    MM    `7                         MM                  
//    MM   d  .gP"Ya `7MMpMMMb.pMMMb.mmMMmm ,pW"Wq.      
//    MM""MM ,M'   Yb  MM    MM    MM  MM  6W'   `Wb   
//    MM   Y 8M""""""  MM    MM    MM  MM  8M     M8 
//    MM     YM.    ,  MM    MM    MM  MM  YA.   ,A9       
//  .JMML.    `Mbmmd'.JMML  JMML  JMML.`Mbmo`Ybmd9'        

void FC_FUNC(g_if_sigma_ccvv_ccov_no0_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO0_X0_TYPE0_NOERI)
  (const double * const W0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no0_x1_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO0_X1_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no1_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO1_X0_TYPE0_NOERI)
  (const double * const W1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no1_x1_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO1_X1_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W1, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no2_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO2_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no2_x1_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO2_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no3_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO3_X0_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const W3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no3_x1_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO3_X1_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W3, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no4_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO4_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no5_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO5_X0_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no6_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO6_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no7_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO7_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no0_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO0_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no0_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO0_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W4, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no1_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO1_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no1_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO1_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W5, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no2_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO2_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no2_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO2_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W6, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no3_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO3_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W7, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no3_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO3_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W7, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no4_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO4_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W8, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no4_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO4_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W8, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no5_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO5_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W9, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no5_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO5_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W9, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no6_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO6_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W10, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no6_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO6_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W10, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no7_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO7_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W11, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no7_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO7_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W11, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no8_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO8_X0_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no9_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO9_X0_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no10_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO10_X0_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no11_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO11_X0_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no12_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO12_X0_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no13_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO13_X0_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no14_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO14_X0_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no15_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO15_X0_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no0_x0_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO0_X0_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const V2, const double * const W24, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no0_x1_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO0_X1_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W24, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no1_x0_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO1_X0_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const V2, const double * const W25, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no1_x1_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO1_X1_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W25, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no2_x0_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO2_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const V2, const double * const W30, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no2_x1_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO2_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W30, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no3_x0_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO3_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const V2, const double * const W31, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no3_x1_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO3_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W31, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no0_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO0_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W12, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no0_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO0_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W12, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no1_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO1_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W13, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no1_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO1_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W13, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no2_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO2_X0_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W14, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no2_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO2_X1_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W14, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no3_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO3_X0_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W15, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no3_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO3_X1_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W15, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no4_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO4_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W16, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no4_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO4_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W16, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no5_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO5_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W17, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no5_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO5_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W17, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no6_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO6_X0_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W18, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no6_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO6_X1_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W18, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no7_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO7_X0_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W19, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no7_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO7_X1_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W19, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no8_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO8_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W20, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no8_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO8_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W20, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no9_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO9_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W21, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no9_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO9_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W21, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no10_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO10_X0_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W22, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no10_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO10_X1_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W22, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no11_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO11_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W23, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no11_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO11_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W23, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no12_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO12_X0_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no13_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO13_X0_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no14_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO14_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no15_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO15_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no16_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO16_X0_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W26, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no16_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO16_X1_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W26, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no17_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO17_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W27, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no17_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO17_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W27, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no18_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO18_X0_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no19_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO19_X0_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no20_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO20_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no21_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO21_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no22_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO22_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W28, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no22_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO22_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W28, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no23_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO23_X0_TYPE1_ERI_V)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W29, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no23_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO23_X1_TYPE1_ERI_V)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W29, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no24_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO24_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W32, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no24_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO24_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W32, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no25_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO25_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W33, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no25_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO25_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W33, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no26_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO26_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W34, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no26_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO26_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W34, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no27_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO27_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W35, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no27_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO27_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W35, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no28_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO28_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no29_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO29_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no30_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO30_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccov_no31_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO31_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 