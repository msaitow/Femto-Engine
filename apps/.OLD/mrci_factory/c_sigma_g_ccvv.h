

#ifndef C_SIGMA_G_CCVV_H
#define C_SIGMA_G_CCVV_H

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

void FC_FUNC(g_if_sigma_g_ccvv_no0_x0, G_IF_SIGMA_G_CCVV_NO0_X0)
  (const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const S0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_g_ccvv_no0_x1, G_IF_SIGMA_G_CCVV_NO0_X1)
  (const FC_INT &sv2, const FC_INT &iv2, 
   const double * const T2, const double * const V2, const double * const S0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

      
 }     
       
       
 #endif
       
       
 