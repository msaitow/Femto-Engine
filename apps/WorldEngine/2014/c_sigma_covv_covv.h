

#ifndef C_SIGMA_COVV_COVV_H
#define C_SIGMA_COVV_COVV_H

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

void FC_FUNC(g_if_sigma_covv_covv_no0_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO0_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const Fc0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no1_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO1_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const Fc0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no2_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO2_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no2_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO2_X1_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no3_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO3_X0_TYPE0_NOERI)
  (const double * const W1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no3_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO3_X1_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W1, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no4_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO4_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no4_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO4_X1_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no5_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO5_X0_TYPE0_NOERI)
  (const double * const W3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no5_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO5_X1_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W3, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no6_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO6_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no6_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO6_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const W4, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no7_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO7_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no7_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO7_X1_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W5, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no8_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO8_X0_TYPE0_NOERI)
  (const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no8_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO8_X1_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const W6, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no9_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO9_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W7, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no9_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO9_X1_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W7, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no0_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO0_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W8, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no0_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO0_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W8, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no1_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO1_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W9, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no1_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO1_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W9, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no2_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO2_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W11, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no2_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO2_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W11, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no3_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO3_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W12, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no3_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO3_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W12, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no4_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO4_X0_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W14, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no4_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO4_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const W14, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no5_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO5_X0_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W20, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no5_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO5_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const W20, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no6_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO6_X0_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W22, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no6_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO6_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const W22, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no7_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO7_X0_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W28, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no7_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO7_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const W28, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no0_x0_type1_eri_o,G_IF_SIGMA_COVV_COVV_NO0_X0_TYPE1_ERI_O)
  (const FC_INT &sa4, const FC_INT &ia4, 
   const double * const V2, const double * const W10, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no1_x0_type1_eri_o,G_IF_SIGMA_COVV_COVV_NO1_X0_TYPE1_ERI_O)
  (const FC_INT &sa4, const FC_INT &ia4, 
   const double * const V2, const double * const W13, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no0_x0_type2_eri_o,G_IF_SIGMA_COVV_COVV_NO0_X0_TYPE2_ERI_O)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W10, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no1_x0_type2_eri_o,G_IF_SIGMA_COVV_COVV_NO1_X0_TYPE2_ERI_O)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W13, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no0_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO0_X0_TYPE1_ERI_V)
  (const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W15, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no0_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO0_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W15, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no1_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO1_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const W16, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no1_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO1_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const W16, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no2_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO2_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W17, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no2_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO2_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W17, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no3_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO3_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const W18, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no3_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO3_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const W18, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no4_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO4_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W19, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no4_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO4_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W19, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no5_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO5_X0_TYPE1_ERI_V)
  (const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W21, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no5_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO5_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W21, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no6_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO6_X0_TYPE1_ERI_V)
  (const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W23, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no6_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO6_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W23, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no7_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO7_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W24, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no7_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO7_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W24, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no8_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO8_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W25, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no8_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO8_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W25, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no9_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO9_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W26, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no9_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO9_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W26, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no10_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO10_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W27, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no10_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO10_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W27, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no11_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO11_X0_TYPE1_ERI_V)
  (const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W29, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no11_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO11_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W29, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no12_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO12_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const W30, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no12_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO12_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W30, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no13_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO13_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W31, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_covv_no13_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO13_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W31, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 