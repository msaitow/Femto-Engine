

#ifndef C_SIGMA_CCOV_COOO_H
#define C_SIGMA_CCOV_COOO_H

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

void FC_FUNC(g_if_sigma_ccov_cooo_no0_x0_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO0_X0_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no0_x1_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO0_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no1_x0_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO1_X0_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no1_x1_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO1_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W1, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no2_x0_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO2_X0_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no2_x1_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO2_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no3_x0_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO3_X0_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no3_x1_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO3_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W3, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no4_x0_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO4_X0_TYPE0_NOERI)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no4_x1_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO4_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const W4, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no5_x0_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO5_X0_TYPE0_NOERI)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no5_x1_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO5_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const W5, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no0_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO0_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W28, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no0_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO0_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W28, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no1_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO1_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W29, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no1_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO1_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W29, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no2_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO2_X0_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, 
   const double * const T2, const double * const W30, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no2_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO2_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa3, const FC_INT &ia3, 
   const double * const V2, const double * const W30, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no3_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO3_X0_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, 
   const double * const T2, const double * const W31, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no3_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO3_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa3, const FC_INT &ia3, 
   const double * const V2, const double * const W31, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no4_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO4_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W32, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no4_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO4_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W32, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no5_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO5_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W33, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no5_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO5_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W33, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no6_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO6_X0_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W36, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no6_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO6_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W36, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no7_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO7_X0_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W37, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no7_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO7_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W37, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no8_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO8_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W38, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no8_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO8_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W38, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no9_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO9_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W39, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no9_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO9_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W39, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no10_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO10_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W40, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no10_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO10_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W40, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no11_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO11_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W41, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no11_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO11_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W41, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no12_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO12_X0_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, 
   const double * const T2, const double * const W42, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no12_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO12_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa3, const FC_INT &ia3, 
   const double * const V2, const double * const W42, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no13_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO13_X0_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, 
   const double * const T2, const double * const W43, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no13_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO13_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa3, const FC_INT &ia3, 
   const double * const V2, const double * const W43, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no14_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO14_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W44, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no14_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO14_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W44, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no15_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO15_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W45, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no15_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO15_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W45, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no16_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO16_X0_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W48, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no16_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO16_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W48, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no17_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO17_X0_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W49, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no17_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO17_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W49, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no18_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO18_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W50, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no18_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO18_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W50, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no19_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO19_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W51, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no19_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO19_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W51, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no0_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO0_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no1_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO1_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W7, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no2_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO2_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W8, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no3_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO3_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W9, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no4_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO4_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W10, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no5_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO5_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W11, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no6_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO6_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W12, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no7_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO7_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W13, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no8_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO8_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W14, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no9_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO9_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W15, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no10_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO10_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W16, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no11_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO11_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W17, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no12_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO12_X0_TYPE0_ERI_V)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W22, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no13_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO13_X0_TYPE0_ERI_V)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W23, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no0_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO0_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W6, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no1_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO1_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W7, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no2_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO2_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W8, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no3_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO3_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W9, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no4_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO4_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W10, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no5_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO5_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W11, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no6_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO6_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W12, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no7_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO7_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W13, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no8_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO8_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W14, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no9_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO9_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W15, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no10_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO10_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W16, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no11_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO11_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W17, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no12_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO12_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const V2, const double * const W18, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no12_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO12_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W18, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no13_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO13_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const V2, const double * const W19, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no13_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO13_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W19, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no14_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO14_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const V2, const double * const W20, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no14_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO14_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W20, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no15_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO15_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const V2, const double * const W21, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no15_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO15_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W21, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no16_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO16_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W22, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no17_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO17_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W23, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no18_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO18_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W24, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no18_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO18_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W24, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no19_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO19_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W25, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no19_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO19_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W25, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no20_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO20_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W26, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no20_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO20_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W26, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no21_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO21_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W27, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no21_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO21_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W27, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no22_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO22_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W34, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no22_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO22_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W34, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no23_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO23_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W35, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no23_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO23_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W35, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no24_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO24_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa1, const FC_INT &ia1, 
   const double * const V2, const double * const W46, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no24_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO24_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W46, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no25_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO25_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa1, const FC_INT &ia1, 
   const double * const V2, const double * const W47, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_cooo_no25_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO25_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W47, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 