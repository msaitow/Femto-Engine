

#ifndef C_SIGMA_OOVV_CCVV_H
#define C_SIGMA_OOVV_CCVV_H

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

void FC_FUNC(g_if_sigma_oovv_ccvv_no0_x0_type1_eri_c,G_IF_SIGMA_OOVV_CCVV_NO0_X0_TYPE1_ERI_C)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const W0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_ccvv_no0_x1_type1_eri_c,G_IF_SIGMA_OOVV_CCVV_NO0_X1_TYPE1_ERI_C)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sc1, const FC_INT &ic1, 
   const double * const T2, const double * const W0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 