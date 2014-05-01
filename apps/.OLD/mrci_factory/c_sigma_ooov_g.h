

#ifndef C_SIGMA_OOOV_G_H
#define C_SIGMA_OOOV_G_H

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

void FC_FUNC(g_if_sigma_ooov_g_no0_x0, G_IF_SIGMA_OOOV_G_NO0_X0)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T0, 
   const double * const h, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_g_y1, G_IF_SIGMA_OOOV_G_Y1)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_g_no0_x1, G_IF_SIGMA_OOOV_G_NO0_X1)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T0, 
   const double * const Y0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_g_y2, G_IF_SIGMA_OOOV_G_Y2)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_g_no0_x2, G_IF_SIGMA_OOOV_G_NO0_X2)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T0, 
   const double * const Y1, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_g_no0_x3, G_IF_SIGMA_OOOV_G_NO0_X3)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_g_no0_x4, G_IF_SIGMA_OOOV_G_NO0_X4)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T0, 
   const double * const V2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

      
 }     
       
       
 #endif
       
       
 