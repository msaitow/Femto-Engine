

#ifndef C_SIGMA_CCOV_CCOV_H
#define C_SIGMA_CCOV_CCOV_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//      _/_/_/_/                            _/             
//     _/        _/_/    _/_/_/  _/_/    _/_/_/_/    _/_/  
//    _/_/_/  _/_/_/_/  _/    _/    _/    _/      _/    _/ 
//   _/      _/        _/    _/    _/    _/      _/    _/  
//  _/        _/_/_/  _/    _/    _/      _/_/    _/_/     

void FC_FUNC(g_if_sigma_ccov_ccov_no0_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO0_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Fc0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no1_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO1_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Fc0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no2_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO2_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Fc0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no3_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO3_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Fc0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no4_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO4_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no4_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO4_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no5_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO5_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no5_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO5_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W1, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no6_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO6_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no6_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO6_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no7_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO7_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no7_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO7_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W3, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no8_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO8_X0_TYPE0_NOERI)
  (const double * const W4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no8_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO8_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W4, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no9_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO9_X0_TYPE0_NOERI)
  (const double * const W5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no9_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO9_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W5, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no10_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO10_X0_TYPE0_NOERI)
  (const double * const W6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no10_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO10_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W6, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no11_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO11_X0_TYPE0_NOERI)
  (const double * const W7, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no11_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO11_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W7, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no12_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO12_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no13_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO13_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no14_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO14_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no15_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO15_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no16_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO16_X0_TYPE0_NOERI)
  (const double * const W8, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no16_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO16_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W8, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no17_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO17_X0_TYPE0_NOERI)
  (const double * const W9, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no17_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO17_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W9, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no18_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO18_X0_TYPE0_NOERI)
  (const double * const W10, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no18_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO18_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W10, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no19_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO19_X0_TYPE0_NOERI)
  (const double * const W11, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no19_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO19_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W11, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no20_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO20_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no21_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO21_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no22_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO22_X0_TYPE0_NOERI)
  (const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W12, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no22_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO22_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const W12, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no23_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO23_X0_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W13, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no23_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO23_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W13, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no24_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO24_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no25_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO25_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no0_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO0_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const V2, const double * const W14, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no0_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO0_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const W14, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no1_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO1_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const V2, const double * const W15, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no1_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO1_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const W15, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no2_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO2_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W16, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no2_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO2_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W16, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no3_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO3_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W17, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no3_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO3_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W17, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no4_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO4_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W18, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no4_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO4_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W18, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no5_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO5_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W19, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no5_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO5_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W19, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no6_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO6_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W20, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no6_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO6_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W20, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no7_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO7_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W21, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no7_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO7_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W21, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no8_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO8_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W22, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no8_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO8_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W22, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no9_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO9_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W23, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no9_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO9_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W23, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no10_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO10_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W26, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no10_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO10_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W26, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no11_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO11_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W27, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no11_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO11_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W27, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no12_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO12_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W28, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no12_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO12_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W28, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no13_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO13_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W29, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no13_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO13_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W29, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no14_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO14_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W30, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no14_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO14_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W30, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no15_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO15_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W31, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no15_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO15_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W31, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no16_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO16_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W32, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no16_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO16_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W32, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no17_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO17_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W33, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no17_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO17_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W33, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no18_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO18_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no19_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO19_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no20_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO20_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W36, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no20_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO20_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W36, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no21_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO21_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W37, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no21_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO21_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W37, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no22_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO22_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W38, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no22_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO22_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W38, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no23_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO23_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W39, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no23_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO23_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W39, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no24_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO24_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W40, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no24_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO24_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W40, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no25_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO25_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W41, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no25_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO25_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W41, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no26_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO26_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W42, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no26_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO26_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W42, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no27_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO27_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W43, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no27_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO27_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W43, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no28_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO28_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W46, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no28_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO28_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W46, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no29_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO29_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W47, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no29_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO29_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W47, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no30_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO30_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W48, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no30_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO30_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W48, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no31_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO31_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W49, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no31_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO31_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W49, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no32_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO32_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W50, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no32_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO32_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W50, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no33_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO33_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W51, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no33_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO33_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W51, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no34_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO34_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W52, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no34_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO34_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W52, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no35_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO35_X0_TYPE1_ERI_C)
  (const FC_INT &sx, const FC_INT &ix, 
   const double * const V2, const double * const W53, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no35_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO35_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const W53, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no36_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO36_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no37_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO37_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no38_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO38_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no39_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO39_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no40_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO40_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no41_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO41_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no42_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO42_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sx, const FC_INT &ix, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no43_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO43_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no0_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO0_X0_TYPE1_ERI_O)
  (const FC_INT &sa4, const FC_INT &ia4, 
   const double * const V2, const double * const W24, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no1_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO1_X0_TYPE1_ERI_O)
  (const FC_INT &sa4, const FC_INT &ia4, 
   const double * const V2, const double * const W25, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no2_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO2_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const V2, const double * const W34, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no2_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO2_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W34, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no3_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO3_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const V2, const double * const W35, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no3_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO3_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W35, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no4_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO4_X0_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, 
   const double * const V2, const double * const W44, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no5_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO5_X0_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, 
   const double * const V2, const double * const W45, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no6_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO6_X0_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const V2, const double * const W54, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no6_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO6_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W54, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no7_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO7_X0_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const V2, const double * const W55, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no7_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO7_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W55, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no8_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO8_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const V2, const double * const W56, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no8_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO8_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W56, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no9_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO9_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const V2, const double * const W57, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no9_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO9_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W57, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no10_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO10_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const V2, const double * const W58, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no10_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO10_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W58, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no11_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO11_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const V2, const double * const W59, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no11_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO11_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W59, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no0_x0_type2_eri_o,G_IF_SIGMA_CCOV_CCOV_NO0_X0_TYPE2_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W24, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no1_x0_type2_eri_o,G_IF_SIGMA_CCOV_CCOV_NO1_X0_TYPE2_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W25, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no2_x0_type2_eri_o,G_IF_SIGMA_CCOV_CCOV_NO2_X0_TYPE2_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W44, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no3_x0_type2_eri_o,G_IF_SIGMA_CCOV_CCOV_NO3_X0_TYPE2_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W45, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no0_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO0_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const V2, const double * const W60, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no0_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO0_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const W60, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no1_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO1_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const V2, const double * const W61, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no1_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO1_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const W61, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no2_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO2_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W62, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no2_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO2_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W62, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no3_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO3_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W63, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no3_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO3_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W63, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no4_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO4_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W64, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no4_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO4_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W64, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no5_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO5_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W65, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no5_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO5_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W65, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no6_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO6_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const V2, const double * const W66, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no6_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO6_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const W66, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no7_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO7_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const V2, const double * const W67, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no7_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO7_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const W67, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no8_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO8_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W68, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no8_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO8_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W68, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no9_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO9_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W69, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no9_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO9_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W69, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no10_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO10_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W70, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no10_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO10_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W70, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no11_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO11_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W71, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no11_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO11_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W71, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no12_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO12_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W72, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no12_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO12_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W72, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no13_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO13_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W73, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no13_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO13_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W73, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no14_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO14_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no15_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO15_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no16_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO16_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no17_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO17_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no18_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO18_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W74, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no18_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO18_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W74, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no19_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO19_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W75, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no19_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO19_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W75, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no20_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO20_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no21_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO21_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no22_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO22_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no23_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO23_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no24_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO24_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W76, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no24_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO24_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W76, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no25_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO25_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W77, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no25_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO25_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W77, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no26_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO26_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W78, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no26_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO26_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W78, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no27_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO27_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W79, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no27_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO27_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W79, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no28_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO28_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W80, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no28_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO28_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W80, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no29_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO29_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W81, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no29_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO29_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W81, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no30_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO30_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W82, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no30_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO30_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W82, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no31_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO31_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W83, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no31_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO31_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W83, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no32_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO32_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no33_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO33_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no34_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO34_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_ccov_no35_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO35_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 