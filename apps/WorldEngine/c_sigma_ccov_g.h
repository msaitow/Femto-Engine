

#ifndef C_SIGMA_CCOV_G_H
#define C_SIGMA_CCOV_G_H

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

void FC_FUNC(g_if_sigma_ccov_g_no0_x0_type1_eri_v,G_IF_SIGMA_CCOV_G_NO0_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_g_no1_x0_type1_eri_v,G_IF_SIGMA_CCOV_G_NO1_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_g_no2_x0_type1_eri_v,G_IF_SIGMA_CCOV_G_NO2_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccov_g_no3_x0_type1_eri_v,G_IF_SIGMA_CCOV_G_NO3_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 