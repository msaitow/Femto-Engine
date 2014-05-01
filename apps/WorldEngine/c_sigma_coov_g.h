

#ifndef C_SIGMA_COOV_G_H
#define C_SIGMA_COOV_G_H

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

void FC_FUNC(g_if_sigma_coov_g_no0_x0_type0_noeri,G_IF_SIGMA_COOV_G_NO0_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T0, 
   const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_g_no0_x0_type1_eri_v,G_IF_SIGMA_COOV_G_NO0_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_g_no1_x0_type1_eri_v,G_IF_SIGMA_COOV_G_NO1_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_g_no2_x0_type1_eri_v,G_IF_SIGMA_COOV_G_NO2_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_coov_g_no3_x0_type1_eri_v,G_IF_SIGMA_COOV_G_NO3_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 