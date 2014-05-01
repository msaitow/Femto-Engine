

#ifndef C_SIGMA_G_OOVV_H
#define C_SIGMA_G_OOVV_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//                                                              
//   _______________                                  ______    
//  |          |                 .'. .`. `````|`````.~      ~.  
//  |______    |______         .'   `   `.    |    |          | 
//  |          |             .'           `.  |    |          | 
//  |          |___________.'               `.|     `.______.'  
//                                                              

void FC_FUNC(g_if_sigma_g_oovv_no0_x0, G_IF_SIGMA_G_OOVV_NO0_X0)
  (const FC_INT &sv2, const FC_INT &iv2, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_g_oovv_no1_x0, G_IF_SIGMA_G_OOVV_NO1_X0)
  (const double * const X, const double * const S0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

      
 }     
       
       
 #endif
       
       
 