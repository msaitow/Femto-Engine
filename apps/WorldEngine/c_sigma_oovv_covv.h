

#ifndef C_SIGMA_OOVV_COVV_H
#define C_SIGMA_OOVV_COVV_H

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

void FC_FUNC(g_if_sigma_oovv_covv_no0_x0_type0_noeri,G_IF_SIGMA_OOVV_COVV_NO0_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no0_x1_type0_noeri,G_IF_SIGMA_OOVV_COVV_NO0_X1_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no1_x0_type0_noeri,G_IF_SIGMA_OOVV_COVV_NO1_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no1_x1_type0_noeri,G_IF_SIGMA_OOVV_COVV_NO1_X1_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W1, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no0_x0_type1_eri_c,G_IF_SIGMA_OOVV_COVV_NO0_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no0_x1_type1_eri_c,G_IF_SIGMA_OOVV_COVV_NO0_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no1_x0_type1_eri_c,G_IF_SIGMA_OOVV_COVV_NO1_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no1_x1_type1_eri_c,G_IF_SIGMA_OOVV_COVV_NO1_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W3, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no0_x0_type1_eri_o,G_IF_SIGMA_OOVV_COVV_NO0_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const W4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no0_x1_type1_eri_o,G_IF_SIGMA_OOVV_COVV_NO0_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sb, const FC_INT &ib, 
   const double * const W4, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no1_x0_type1_eri_o,G_IF_SIGMA_OOVV_COVV_NO1_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const W8, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no1_x1_type1_eri_o,G_IF_SIGMA_OOVV_COVV_NO1_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sb, const FC_INT &ib, 
   const double * const W8, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no0_x0_type1_eri_v,G_IF_SIGMA_OOVV_COVV_NO0_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const W5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no0_x1_type1_eri_v,G_IF_SIGMA_OOVV_COVV_NO0_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const W5, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no1_x0_type1_eri_v,G_IF_SIGMA_OOVV_COVV_NO1_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const W6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no1_x1_type1_eri_v,G_IF_SIGMA_OOVV_COVV_NO1_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const W6, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no2_x0_type1_eri_v,G_IF_SIGMA_OOVV_COVV_NO2_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const W7, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no2_x1_type1_eri_v,G_IF_SIGMA_OOVV_COVV_NO2_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const W7, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no3_x0_type1_eri_v,G_IF_SIGMA_OOVV_COVV_NO3_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W9, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no3_x1_type1_eri_v,G_IF_SIGMA_OOVV_COVV_NO3_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W9, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no4_x0_type1_eri_v,G_IF_SIGMA_OOVV_COVV_NO4_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W10, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no4_x1_type1_eri_v,G_IF_SIGMA_OOVV_COVV_NO4_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W10, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no5_x0_type1_eri_v,G_IF_SIGMA_OOVV_COVV_NO5_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const W11, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_covv_no5_x1_type1_eri_v,G_IF_SIGMA_OOVV_COVV_NO5_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sb, const FC_INT &ib, 
   const double * const W11, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 