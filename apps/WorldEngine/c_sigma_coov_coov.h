

#ifndef C_SIGMA_COOV_COOV_H
#define C_SIGMA_COOV_COOV_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//  `MMMMMMM                                         
//   MM    \                         /              
//   MM       ____  ___  __    __   /M      _____    
//   MM   ,  6MMMMb `MM 6MMb  6MMb /MMMMM  6MMMMMb   
//   MMMMMM 6M'  `Mb MM69 `MM69 `Mb MM    6M'   `Mb  
//   MM   ` MM    MM MM'   MM'   MM MM    MM     MM  
//   MM     MMMMMMMM MM    MM    MM MM    MM     MM  
//   MM     MM       MM    MM    MM MM    MM     MM  
//   MM     YM    d9 MM    MM    MM YM.  ,YM.   ,M9  
//  _MM_     YMMMM9 _MM_  _MM_  _MM_ YMMM9 YMMMMM9   

void FC_FUNC(g_if_sigma_coov_coov_no0_x0_type0_noeri,G_IF_SIGMA_COOV_COOV_NO0_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Fc0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no1_x0_type0_noeri,G_IF_SIGMA_COOV_COOV_NO1_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const Fc0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no2_x0_type0_noeri,G_IF_SIGMA_COOV_COOV_NO2_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no2_x1_type0_noeri,G_IF_SIGMA_COOV_COOV_NO2_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no3_x0_type0_noeri,G_IF_SIGMA_COOV_COOV_NO3_X0_TYPE0_NOERI)
  (const double * const W1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no3_x1_type0_noeri,G_IF_SIGMA_COOV_COOV_NO3_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W1, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no4_x0_type0_noeri,G_IF_SIGMA_COOV_COOV_NO4_X0_TYPE0_NOERI)
  (const double * const W2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no4_x1_type0_noeri,G_IF_SIGMA_COOV_COOV_NO4_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no5_x0_type0_noeri,G_IF_SIGMA_COOV_COOV_NO5_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no5_x1_type0_noeri,G_IF_SIGMA_COOV_COOV_NO5_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W3, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no6_x0_type0_noeri,G_IF_SIGMA_COOV_COOV_NO6_X0_TYPE0_NOERI)
  (const double * const W4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no6_x1_type0_noeri,G_IF_SIGMA_COOV_COOV_NO6_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W4, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no7_x0_type0_noeri,G_IF_SIGMA_COOV_COOV_NO7_X0_TYPE0_NOERI)
  (const double * const W5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no7_x1_type0_noeri,G_IF_SIGMA_COOV_COOV_NO7_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W5, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no8_x0_type0_noeri,G_IF_SIGMA_COOV_COOV_NO8_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no8_x1_type0_noeri,G_IF_SIGMA_COOV_COOV_NO8_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W6, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no9_x0_type0_noeri,G_IF_SIGMA_COOV_COOV_NO9_X0_TYPE0_NOERI)
  (const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W7, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no9_x1_type0_noeri,G_IF_SIGMA_COOV_COOV_NO9_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const W7, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no10_x0_type0_noeri,G_IF_SIGMA_COOV_COOV_NO10_X0_TYPE0_NOERI)
  (const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W8, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no10_x1_type0_noeri,G_IF_SIGMA_COOV_COOV_NO10_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const W8, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no11_x0_type0_noeri,G_IF_SIGMA_COOV_COOV_NO11_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no0_x0_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO0_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W9, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no0_x1_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO0_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W9, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no1_x0_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO1_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W10, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no1_x1_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO1_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W10, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no2_x0_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO2_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W11, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no2_x1_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO2_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W11, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no3_x0_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO3_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W12, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no3_x1_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO3_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W12, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no4_x0_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO4_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W14, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no4_x1_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO4_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W14, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no5_x0_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO5_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W15, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no5_x1_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO5_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W15, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no6_x0_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO6_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W17, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no6_x1_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO6_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const W17, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no7_x0_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO7_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W18, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no7_x1_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO7_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const W18, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no8_x0_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO8_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W20, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no8_x1_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO8_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const W20, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no9_x0_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO9_X0_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W22, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no9_x1_type1_eri_c,G_IF_SIGMA_COOV_COOV_NO9_X1_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sw, const FC_INT &iw, 
   const double * const W22, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no0_x0_type1_eri_o,G_IF_SIGMA_COOV_COOV_NO0_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W13, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no0_x1_type1_eri_o,G_IF_SIGMA_COOV_COOV_NO0_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W13, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no1_x0_type1_eri_o,G_IF_SIGMA_COOV_COOV_NO1_X0_TYPE1_ERI_O)
  (const FC_INT &sa4, const FC_INT &ia4, 
   const double * const V2, const double * const W16, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no2_x0_type1_eri_o,G_IF_SIGMA_COOV_COOV_NO2_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const V2, const double * const W19, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no2_x1_type1_eri_o,G_IF_SIGMA_COOV_COOV_NO2_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W19, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no3_x0_type1_eri_o,G_IF_SIGMA_COOV_COOV_NO3_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W21, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no3_x1_type1_eri_o,G_IF_SIGMA_COOV_COOV_NO3_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W21, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no4_x0_type1_eri_o,G_IF_SIGMA_COOV_COOV_NO4_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W23, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no4_x1_type1_eri_o,G_IF_SIGMA_COOV_COOV_NO4_X1_TYPE1_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W23, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no0_x0_type2_eri_o,G_IF_SIGMA_COOV_COOV_NO0_X0_TYPE2_ERI_O)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W16, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no0_x0_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO0_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W24, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no0_x1_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO0_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W24, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no1_x0_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO1_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W25, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no1_x1_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO1_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W25, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no2_x0_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO2_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W26, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no2_x1_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO2_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W26, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no3_x0_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO3_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W27, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no3_x1_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO3_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W27, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no4_x0_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO4_X0_TYPE1_ERI_V)
  (const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W28, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no4_x1_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO4_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W28, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no5_x0_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO5_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W29, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no5_x1_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO5_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W29, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no6_x0_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO6_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W30, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no6_x1_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO6_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W30, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no7_x0_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO7_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W31, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no7_x1_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO7_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W31, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no8_x0_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO8_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W32, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no8_x1_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO8_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W32, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no9_x0_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO9_X0_TYPE1_ERI_V)
  (const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W33, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no9_x1_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO9_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W33, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no10_x0_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO10_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W34, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no10_x1_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO10_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W34, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no11_x0_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO11_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W35, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no11_x1_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO11_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W35, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no12_x0_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO12_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W36, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no12_x1_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO12_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W36, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no13_x0_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO13_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W37, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_coov_no13_x1_type1_eri_v,G_IF_SIGMA_COOV_COOV_NO13_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W37, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 