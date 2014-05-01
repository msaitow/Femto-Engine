

#ifndef C_SIGMA_CCVV_CCOO_H
#define C_SIGMA_CCVV_CCOO_H

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

void FC_FUNC(g_if_sigma_ccvv_ccoo_no0_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOO_NO0_X0_TYPE1_ERI_V)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccoo_no0_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOO_NO0_X1_TYPE1_ERI_V)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccoo_no1_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOO_NO1_X0_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccoo_no1_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOO_NO1_X1_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W1, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccoo_no2_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOO_NO2_X0_TYPE1_ERI_V)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccoo_no2_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOO_NO2_X1_TYPE1_ERI_V)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccoo_no3_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOO_NO3_X0_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccoo_no3_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOO_NO3_X1_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W3, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccoo_no4_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOO_NO4_X0_TYPE1_ERI_V)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const W4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccoo_no4_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOO_NO4_X1_TYPE1_ERI_V)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W4, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccoo_no5_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOO_NO5_X0_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ccvv_ccoo_no6_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOO_NO6_X0_TYPE1_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 