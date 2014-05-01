

#ifndef C_SIGMA_G_CCOO_H
#define C_SIGMA_G_CCOO_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//  8888888888                     888                  
//  888                            888                  
//  888                            888                  
//  8888888  .d88b.  88888b.d88b.  888888  .d88b.       
//  888     d8P  Y8b 888 "888 "88b 888    d88""88b  
//  888     88888888 888  888  888 888    888  888      
//  888     Y8b.     888  888  888 Y88b.  Y88..88P      
//  888      "Y8888  888  888  888  "Y888  "Y88P"   

void FC_FUNC(g_if_sigma_g_ccoo_no0_x0_type1_eri_c,G_IF_SIGMA_G_CCOO_NO0_X0_TYPE1_ERI_C)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const W0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_g_ccoo_no0_x1_type1_eri_c,G_IF_SIGMA_G_CCOO_NO0_X1_TYPE1_ERI_C)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sc1, const FC_INT &ic1, 
   const double * const T2, const double * const W0, const double * const S0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_g_ccoo_no1_x0_type1_eri_c,G_IF_SIGMA_G_CCOO_NO1_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const W1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_g_ccoo_no1_x1_type1_eri_c,G_IF_SIGMA_G_CCOO_NO1_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc1, const FC_INT &ic1, 
   const double * const T2, const double * const W1, const double * const S0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_g_ccoo_no2_x0_type1_eri_c,G_IF_SIGMA_G_CCOO_NO2_X0_TYPE1_ERI_C)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const W2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_g_ccoo_no2_x1_type1_eri_c,G_IF_SIGMA_G_CCOO_NO2_X1_TYPE1_ERI_C)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sc1, const FC_INT &ic1, 
   const double * const T2, const double * const W2, const double * const S0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_g_ccoo_no3_x0_type1_eri_c,G_IF_SIGMA_G_CCOO_NO3_X0_TYPE1_ERI_C)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sc1, const FC_INT &ic1, 
   const double * const T2, const double * const V2, const double * const S0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_g_ccoo_no4_x0_type1_eri_c,G_IF_SIGMA_G_CCOO_NO4_X0_TYPE1_ERI_C)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sc1, const FC_INT &ic1, 
   const double * const T2, const double * const V2, const double * const S0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 