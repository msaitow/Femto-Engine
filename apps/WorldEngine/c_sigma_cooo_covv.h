

#ifndef C_SIGMA_COOO_COVV_H
#define C_SIGMA_COOO_COVV_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//  ___________                __               
//  \_   _____/____    _____ _/  |_  ____      
//   |    __)_/ __ \  /     \\   __\/  _ \ 
//   |     \ \  ___/ |  Y Y  \|  | (  <_> )  
//   \___  /  \___  >|__|_|  /|__|  \____/   
//       \/       \/       \/                

void FC_FUNC(g_if_sigma_cooo_covv_no0_x0_type1_eri_o,G_IF_SIGMA_COOO_COVV_NO0_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_covv_no0_x1_type1_eri_o,G_IF_SIGMA_COOO_COVV_NO0_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W3, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_covv_no1_x0_type1_eri_o,G_IF_SIGMA_COOO_COVV_NO1_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const W4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_covv_no1_x1_type1_eri_o,G_IF_SIGMA_COOO_COVV_NO1_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W4, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_covv_no2_x0_type1_eri_o,G_IF_SIGMA_COOO_COVV_NO2_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const W5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_covv_no2_x1_type1_eri_o,G_IF_SIGMA_COOO_COVV_NO2_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W5, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_covv_no3_x0_type1_eri_o,G_IF_SIGMA_COOO_COVV_NO3_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_covv_no3_x1_type1_eri_o,G_IF_SIGMA_COOO_COVV_NO3_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W6, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_covv_no0_x0_type1_eri_v,G_IF_SIGMA_COOO_COVV_NO0_X0_TYPE1_ERI_V)
  (const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const W0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_covv_no1_x0_type1_eri_v,G_IF_SIGMA_COOO_COVV_NO1_X0_TYPE1_ERI_V)
  (const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const W1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_covv_no2_x0_type1_eri_v,G_IF_SIGMA_COOO_COVV_NO2_X0_TYPE1_ERI_V)
  (const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const W2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_covv_no0_x0_type2_eri_v,G_IF_SIGMA_COOO_COVV_NO0_X0_TYPE2_ERI_V)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_covv_no1_x0_type2_eri_v,G_IF_SIGMA_COOO_COVV_NO1_X0_TYPE2_ERI_V)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W1, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_covv_no2_x0_type2_eri_v,G_IF_SIGMA_COOO_COVV_NO2_X0_TYPE2_ERI_V)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 