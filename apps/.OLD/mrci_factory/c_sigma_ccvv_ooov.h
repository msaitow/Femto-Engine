

#ifndef C_SIGMA_CCVV_OOOV_H
#define C_SIGMA_CCVV_OOOV_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//  `MMMMMMM                                         
//   MM    \                         /              
//   MM       ____  ___  __    __   /M      _____    
//   MM   ,  6MMMMb `MM 6MMb  6MMb /MMMMM  6MMMMMb   
//   MMMMMM 6M'  `Mb MM69 `MM69 `Mb MM    6M'   `Mb  
//   MM   ` MM    MM MM'   MM'   MM MM    MM     MM  
//   MM     MMMMMMMM MM    MM    MM MM    MM     MM  
//   MM     MM       MM    MM    MM MM    MM     MM  
//   MM     YM    d9 MM    MM    MM YM.  ,YM.   ,M9  
//  _MM_     YMMMM9 _MM_  _MM_  _MM_ YMMM9 YMMMMM9   

void FC_FUNC(g_if_sigma_ccvv_ooov_no0_x0, G_IF_SIGMA_CCVV_OOOV_NO0_X0)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ccvv_ooov_no1_x0, G_IF_SIGMA_CCVV_OOOV_NO1_X0)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ccvv_ooov_no0_x1, G_IF_SIGMA_CCVV_OOOV_NO0_X1)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ccvv_ooov_no1_x1, G_IF_SIGMA_CCVV_OOOV_NO1_X1)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ccvv_ooov_no0_x2, G_IF_SIGMA_CCVV_OOOV_NO0_X2)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ccvv_ooov_no1_x2, G_IF_SIGMA_CCVV_OOOV_NO1_X2)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ccvv_ooov_no0_x3, G_IF_SIGMA_CCVV_OOOV_NO0_X3)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ccvv_ooov_no1_x3, G_IF_SIGMA_CCVV_OOOV_NO1_X3)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

      
 }     
       
       
 #endif
       
       
 