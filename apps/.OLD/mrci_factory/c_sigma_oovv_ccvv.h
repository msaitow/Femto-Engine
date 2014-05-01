

#ifndef C_SIGMA_OOVV_CCVV_H
#define C_SIGMA_OOVV_CCVV_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//        :::::::::: ::::::::::   :::   ::: ::::::::::: :::::::: 
//       :+:        :+:         :+:+: :+:+:    :+:    :+:    :+: 
//      +:+        +:+        +:+ +:+:+ +:+   +:+    +:+    +:+  
//     :#::+::#   +#++:++#   +#+  +:+  +#+   +#+    +#+    +:+   
//    +#+        +#+        +#+       +#+   +#+    +#+    +#+    
//   #+#        #+#        #+#       #+#   #+#    #+#    #+#     
//  ###        ########## ###       ###   ###     ########       

void FC_FUNC(g_if_sigma_oovv_ccvv_no0_x0, G_IF_SIGMA_OOVV_CCVV_NO0_X0)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sc1, const FC_INT &ic1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ccvv_no1_x0, G_IF_SIGMA_OOVV_CCVV_NO1_X0)
  (const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ccvv_no0_x1, G_IF_SIGMA_OOVV_CCVV_NO0_X1)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, const FC_INT &sc1, const FC_INT &ic1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_ccvv_no1_x1, G_IF_SIGMA_OOVV_CCVV_NO1_X1)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc, const FC_INT &ic, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

      
 }     
       
       
 #endif
       
       
 