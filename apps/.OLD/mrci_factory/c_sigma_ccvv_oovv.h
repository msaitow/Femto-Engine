

#ifndef C_SIGMA_CCVV_OOVV_H
#define C_SIGMA_CCVV_OOVV_H

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

void FC_FUNC(g_if_sigma_ccvv_oovv_no0_x0, G_IF_SIGMA_CCVV_OOVV_NO0_X0)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ccvv_oovv_no1_x0, G_IF_SIGMA_CCVV_OOVV_NO1_X0)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ccvv_oovv_no0_x1, G_IF_SIGMA_CCVV_OOVV_NO0_X1)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ccvv_oovv_no1_x1, G_IF_SIGMA_CCVV_OOVV_NO1_X1)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

      
 }     
       
       
 #endif
       
       
 