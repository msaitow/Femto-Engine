

#ifndef C_SIGMA_COVV_OOVV_H
#define C_SIGMA_COVV_OOVV_H

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

void FC_FUNC(g_if_sigma_covv_oovv_no0_x0_type0_noeri,G_IF_SIGMA_COVV_OOVV_NO0_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_oovv_no0_x1_type0_noeri,G_IF_SIGMA_COVV_OOVV_NO0_X1_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_oovv_no0_x0_type1_eri_c,G_IF_SIGMA_COVV_OOVV_NO0_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_oovv_no0_x1_type1_eri_c,G_IF_SIGMA_COVV_OOVV_NO0_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W1, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_oovv_no0_x0_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO0_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_oovv_no0_x1_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO0_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_oovv_no1_x0_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO1_X0_TYPE1_ERI_V)
  (const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_oovv_no1_x1_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO1_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W3, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_oovv_no2_x0_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO2_X0_TYPE1_ERI_V)
  (const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_oovv_no2_x1_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO2_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W4, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_oovv_no3_x0_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO3_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_covv_oovv_no3_x1_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO3_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W5, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 