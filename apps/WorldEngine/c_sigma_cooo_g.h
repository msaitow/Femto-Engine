

#ifndef C_SIGMA_COOO_G_H
#define C_SIGMA_COOO_G_H

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

void FC_FUNC(g_if_sigma_cooo_g_no0_x0_type0_noeri,G_IF_SIGMA_COOO_G_NO0_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T0, 
   const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_g_no1_x0_type0_noeri,G_IF_SIGMA_COOO_G_NO1_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T0, 
   const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_g_no2_x0_type0_noeri,G_IF_SIGMA_COOO_G_NO2_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T0, 
   const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_g_no0_x0_type1_eri_c,G_IF_SIGMA_COOO_G_NO0_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_g_no1_x0_type1_eri_c,G_IF_SIGMA_COOO_G_NO1_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_g_no2_x0_type1_eri_c,G_IF_SIGMA_COOO_G_NO2_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_g_no0_x0_type1_eri_o,G_IF_SIGMA_COOO_G_NO0_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_g_no1_x0_type1_eri_o,G_IF_SIGMA_COOO_G_NO1_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_g_no2_x0_type1_eri_o,G_IF_SIGMA_COOO_G_NO2_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_g_no3_x0_type1_eri_o,G_IF_SIGMA_COOO_G_NO3_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 