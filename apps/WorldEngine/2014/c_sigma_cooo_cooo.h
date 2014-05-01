

#ifndef C_SIGMA_COOO_COOO_H
#define C_SIGMA_COOO_COOO_H

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

void FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no0_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO0_X1_TYPE0_NOERI)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no1_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO1_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W1, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE0_NOERI)
  (const double * const W2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no2_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO2_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no3_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO3_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W3, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE0_NOERI)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no4_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO4_X1_TYPE0_NOERI)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const W4, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no5_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO5_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W5, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no6_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO6_X0_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const W6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no6_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO6_X1_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W6, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no7_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO7_X0_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const W7, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no7_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO7_X1_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W7, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no8_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO8_X0_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W8, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no8_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO8_X1_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W8, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no9_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO9_X0_TYPE0_NOERI)
  (const double * const W9, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no9_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO9_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W9, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no10_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO10_X0_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W10, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no10_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO10_X1_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W10, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no11_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO11_X0_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W11, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no11_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO11_X1_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W11, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no12_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO12_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W12, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no12_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO12_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W12, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no13_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO13_X0_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const W13, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no13_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO13_X1_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W13, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no14_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO14_X0_TYPE0_NOERI)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W14, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no14_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO14_X1_TYPE0_NOERI)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const W14, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no15_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO15_X0_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W15, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no15_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO15_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W15, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no16_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO16_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W16, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no16_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO16_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W16, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no17_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO17_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W17, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no17_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO17_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W17, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no18_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO18_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W18, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no18_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO18_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W18, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no19_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO19_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W19, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no19_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO19_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W19, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no20_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO20_X0_TYPE0_NOERI)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W20, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no20_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO20_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W20, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no21_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO21_X0_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const W21, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no21_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO21_X1_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W21, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no22_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO22_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W22, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no22_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO22_X1_TYPE0_NOERI)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W22, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no23_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO23_X0_TYPE0_NOERI)
  (const double * const W23, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no23_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO23_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W23, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no24_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO24_X0_TYPE0_NOERI)
  (const double * const W24, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no24_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO24_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W24, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no25_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO25_X0_TYPE0_NOERI)
  (const double * const W25, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no25_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO25_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W25, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no26_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO26_X0_TYPE0_NOERI)
  (const double * const W26, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no26_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO26_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W26, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no27_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO27_X0_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const W27, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no27_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO27_X1_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W27, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no28_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO28_X0_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const E0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no29_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO29_X0_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const E0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no30_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO30_X0_TYPE0_NOERI)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const E0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no31_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO31_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const E0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no32_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO32_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const E0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no33_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO33_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const E0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no34_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO34_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const E0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no35_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO35_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W403, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no35_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO35_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W403, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no36_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO36_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W404, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no36_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO36_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W404, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no37_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO37_X0_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W405, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no37_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO37_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W405, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no38_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO38_X0_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W406, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no38_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO38_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W406, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no39_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO39_X0_TYPE0_NOERI)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W407, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no39_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO39_X1_TYPE0_NOERI)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const W407, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no40_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO40_X0_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W408, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no40_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO40_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W408, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no41_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO41_X0_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W409, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no41_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO41_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W409, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no42_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO42_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W410, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no42_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO42_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W410, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no43_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO43_X0_TYPE0_NOERI)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W411, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no43_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO43_X1_TYPE0_NOERI)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const W411, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no44_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO44_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W412, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no44_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO44_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W412, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no45_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO45_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W413, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no45_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO45_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W413, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no46_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO46_X0_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W414, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no46_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO46_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W414, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no47_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO47_X0_TYPE0_NOERI)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W415, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no47_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO47_X1_TYPE0_NOERI)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W415, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W135, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W137, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W143, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W147, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W151, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W161, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no6_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO6_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W163, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no7_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO7_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W167, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no8_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO8_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W169, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no9_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO9_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W181, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no10_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO10_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W183, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no11_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO11_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W185, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no12_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO12_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W189, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no13_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO13_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W191, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no14_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO14_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W193, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no15_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO15_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W195, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no16_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO16_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W197, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no17_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO17_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W201, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no18_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO18_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W205, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no19_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO19_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W215, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no20_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO20_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W217, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no21_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO21_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W221, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no22_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO22_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W223, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no23_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO23_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W235, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no24_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO24_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W237, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no25_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO25_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W239, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no26_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO26_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W243, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no27_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO27_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W245, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no28_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO28_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W247, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no29_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO29_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W249, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no30_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO30_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W275, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no31_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO31_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W277, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no32_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO32_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W283, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no33_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO33_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W285, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no34_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO34_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W289, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no35_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO35_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W293, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no36_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO36_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W295, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no37_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO37_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W297, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no38_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO38_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W299, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no39_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO39_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W303, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no40_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO40_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W305, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no41_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO41_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W311, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no42_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO42_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W313, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no43_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO43_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W315, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no44_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO44_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W317, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no45_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO45_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W321, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no46_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO46_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W323, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no47_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO47_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W325, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no48_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO48_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W329, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no49_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO49_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W331, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no50_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO50_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W333, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no51_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO51_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W335, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no52_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO52_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W337, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no53_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO53_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W339, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no54_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO54_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W343, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no55_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO55_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W347, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no56_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO56_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W349, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no57_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO57_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W351, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no58_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO58_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W353, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no59_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO59_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W357, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no60_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO60_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W359, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no61_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO61_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W365, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no62_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO62_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W367, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no63_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO63_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W369, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no64_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO64_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W371, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no65_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO65_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W375, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no66_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO66_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W377, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no67_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO67_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W379, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no68_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO68_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W383, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no69_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO69_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W385, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no70_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO70_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W387, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no71_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO71_X0_TYPE0_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W389, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W37, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no0_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO0_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W37, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W38, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no1_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO1_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W38, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W42, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no2_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO2_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W42, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W50, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no3_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO3_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W50, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W63, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no4_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO4_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W63, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W64, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no5_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO5_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W64, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no6_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO6_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W68, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no6_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO6_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W68, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no7_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO7_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W76, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no7_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO7_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W76, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no8_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO8_X0_TYPE1_ERI_C)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W87, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no8_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO8_X1_TYPE1_ERI_C)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W87, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no9_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO9_X0_TYPE1_ERI_C)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W88, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no9_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO9_X1_TYPE1_ERI_C)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W88, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no10_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO10_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W89, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no10_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO10_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W89, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no11_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO11_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W90, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no11_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO11_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W90, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no12_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO12_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W91, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no12_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO12_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W91, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no13_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO13_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W92, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no13_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO13_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W92, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no14_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO14_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W93, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no14_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO14_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W93, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no15_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO15_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W94, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no15_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO15_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W94, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no16_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO16_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W97, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no16_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO16_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W97, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no17_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO17_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W98, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no17_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO17_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W98, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no18_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO18_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W101, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no18_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO18_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W101, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no19_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO19_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W102, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no19_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO19_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W102, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no20_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO20_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W103, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no20_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO20_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W103, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no21_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO21_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W104, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no21_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO21_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W104, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no22_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO22_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W105, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no22_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO22_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W105, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no23_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO23_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W106, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no23_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO23_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W106, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no24_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO24_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W107, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no24_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO24_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W107, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no25_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO25_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W108, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no25_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO25_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W108, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no26_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO26_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W113, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no26_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO26_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W113, const double * const W114, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no26_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO26_X2_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W114, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no27_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO27_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W115, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no27_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO27_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W115, const double * const W116, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no27_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO27_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W116, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no28_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO28_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W118, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no28_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO28_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W117, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no28_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO28_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W117, const double * const W118, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no29_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO29_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W119, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no29_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO29_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W120, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no29_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO29_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W119, const double * const W120, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no30_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO30_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W121, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no30_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO30_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const W121, const double * const W122, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no30_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO30_X2_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W122, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no31_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO31_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W123, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no31_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO31_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const W123, const double * const W124, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no31_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO31_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W124, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no32_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO32_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W127, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no32_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO32_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W127, const double * const W128, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no32_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO32_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W128, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no33_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO33_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W129, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no33_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO33_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W129, const double * const W130, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no33_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO33_X2_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W130, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no34_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO34_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W131, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no34_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO34_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W132, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no34_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO34_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W131, const double * const W132, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no35_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO35_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W135, const double * const W136, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no35_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO35_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W136, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no36_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO36_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W137, const double * const W138, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no36_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO36_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W138, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no37_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO37_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W139, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no37_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO37_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W139, const double * const W140, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no37_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO37_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W140, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no38_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO38_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W141, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no38_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO38_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const W141, const double * const W142, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no38_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO38_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W142, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no39_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO39_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W143, const double * const W144, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no39_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO39_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W144, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no40_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO40_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W147, const double * const W148, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no40_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO40_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W148, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no41_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO41_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W151, const double * const W152, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no41_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO41_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W152, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no42_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO42_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W153, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no42_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO42_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W154, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no42_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO42_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W153, const double * const W154, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no43_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO43_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W155, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no43_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO43_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W156, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no43_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO43_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W155, const double * const W156, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no44_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO44_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W157, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no44_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO44_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W158, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no44_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO44_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W157, const double * const W158, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no45_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO45_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W159, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no45_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO45_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W160, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no45_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO45_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W159, const double * const W160, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no46_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO46_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W161, const double * const W162, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no46_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO46_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W162, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no47_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO47_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W163, const double * const W164, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no47_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO47_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W164, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no48_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO48_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W166, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no48_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO48_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W165, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no48_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO48_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W165, const double * const W166, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no49_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO49_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W167, const double * const W168, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no49_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO49_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W168, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no50_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO50_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W169, const double * const W170, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no50_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO50_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W170, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no51_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO51_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W172, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no51_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO51_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W171, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no51_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO51_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W171, const double * const W172, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no52_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO52_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W174, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no52_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO52_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W173, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no52_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO52_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W173, const double * const W174, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no53_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO53_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W176, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no53_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO53_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W175, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no53_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO53_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W175, const double * const W176, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no54_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO54_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W178, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no54_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO54_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W177, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no54_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO54_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W177, const double * const W178, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no55_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO55_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W179, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no55_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO55_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W180, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no55_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO55_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W179, const double * const W180, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no56_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO56_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W181, const double * const W182, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no56_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO56_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W182, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no57_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO57_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W183, const double * const W184, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no57_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO57_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W184, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no58_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO58_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W185, const double * const W186, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no58_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO58_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W186, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no59_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO59_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W187, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no59_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO59_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W188, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no59_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO59_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W187, const double * const W188, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no60_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO60_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W189, const double * const W190, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no60_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO60_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W190, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no61_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO61_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W191, const double * const W192, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no61_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO61_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W192, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no62_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO62_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W193, const double * const W194, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no62_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO62_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W194, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no63_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO63_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W195, const double * const W196, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no63_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO63_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W196, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no64_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO64_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W197, const double * const W198, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no64_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO64_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W198, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no65_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO65_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W201, const double * const W202, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no65_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO65_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W202, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no66_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO66_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W205, const double * const W206, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no66_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO66_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W206, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no67_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO67_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W207, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no67_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO67_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W208, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no67_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO67_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W207, const double * const W208, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no68_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO68_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W209, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no68_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO68_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W210, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no68_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO68_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W209, const double * const W210, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no69_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO69_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W211, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no69_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO69_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W212, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no69_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO69_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W211, const double * const W212, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no70_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO70_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W213, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no70_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO70_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W214, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no70_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO70_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W213, const double * const W214, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no71_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO71_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W215, const double * const W216, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no71_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO71_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W216, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no72_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO72_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W217, const double * const W218, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no72_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO72_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W218, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no73_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO73_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W220, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no73_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO73_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W219, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no73_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO73_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W219, const double * const W220, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no74_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO74_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W221, const double * const W222, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no74_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO74_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W222, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no75_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO75_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W223, const double * const W224, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no75_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO75_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W224, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no76_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO76_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W226, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no76_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO76_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W225, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no76_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO76_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W225, const double * const W226, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no77_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO77_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W228, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no77_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO77_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W227, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no77_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO77_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W227, const double * const W228, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no78_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO78_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W230, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no78_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO78_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W229, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no78_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO78_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W229, const double * const W230, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no79_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO79_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W232, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no79_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO79_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W231, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no79_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO79_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W231, const double * const W232, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no80_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO80_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W233, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no80_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO80_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W234, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no80_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO80_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W233, const double * const W234, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no81_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO81_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W235, const double * const W236, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no81_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO81_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W236, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no82_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO82_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W237, const double * const W238, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no82_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO82_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W238, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no83_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO83_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W239, const double * const W240, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no83_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO83_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W240, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no84_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO84_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W241, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no84_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO84_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W242, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no84_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO84_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W241, const double * const W242, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no85_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO85_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W243, const double * const W244, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no85_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO85_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W244, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no86_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO86_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W245, const double * const W246, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no86_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO86_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W246, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no87_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO87_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W247, const double * const W248, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no87_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO87_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W248, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no88_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO88_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W249, const double * const W250, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no88_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO88_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W250, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no89_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO89_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W251, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no89_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO89_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W251, const double * const W252, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no90_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO90_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W253, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no90_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO90_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W253, const double * const W254, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no90_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO90_X2_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W254, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no91_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO91_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W255, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no91_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO91_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W255, const double * const W256, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no91_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO91_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W256, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no92_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO92_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W258, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no92_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO92_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W257, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no92_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO92_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W257, const double * const W258, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no93_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO93_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W259, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no93_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO93_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W260, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no93_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO93_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W259, const double * const W260, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no94_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO94_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W261, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no94_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO94_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W261, const double * const W262, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no94_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO94_X2_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W262, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no95_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO95_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W263, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no95_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO95_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W263, const double * const W264, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no95_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO95_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W264, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no96_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO96_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W265, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no96_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO96_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W266, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no96_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO96_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W265, const double * const W266, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no97_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO97_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W267, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no97_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO97_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const W267, const double * const W268, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no97_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO97_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W268, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no98_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO98_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W269, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no98_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO98_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const W269, const double * const W270, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no98_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO98_X2_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W270, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no99_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO99_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W271, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no99_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO99_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W272, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no99_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO99_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W271, const double * const W272, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no100_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO100_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W273, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no100_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO100_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W273, const double * const W274, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no101_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO101_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W275, const double * const W276, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no101_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO101_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W276, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no102_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO102_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W277, const double * const W278, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no102_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO102_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W278, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no103_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO103_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W279, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no103_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO103_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const W279, const double * const W280, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no103_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO103_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W280, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no104_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO104_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const T2, const double * const V2, const double * const W281, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no104_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO104_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sw, const FC_INT &iw, 
   const double * const W281, const double * const W282, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no104_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO104_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W282, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no105_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO105_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W283, const double * const W284, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no105_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO105_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W284, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no106_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO106_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W285, const double * const W286, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no106_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO106_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W286, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no107_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO107_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W288, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no107_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO107_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W287, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no107_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO107_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W287, const double * const W288, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no108_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO108_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W289, const double * const W290, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no108_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO108_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W290, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no109_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO109_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W292, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no109_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO109_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W291, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no109_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO109_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W291, const double * const W292, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no110_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO110_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W293, const double * const W294, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no110_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO110_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W294, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no111_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO111_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W295, const double * const W296, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no111_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO111_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W296, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no112_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO112_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W297, const double * const W298, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no112_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO112_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W298, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no113_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO113_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W299, const double * const W300, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no113_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO113_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W300, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no114_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO114_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W301, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no114_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO114_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W301, const double * const W302, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no115_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO115_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W303, const double * const W304, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no115_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO115_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W304, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no116_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO116_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W305, const double * const W306, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no116_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO116_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W306, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no117_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO117_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W307, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no117_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO117_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W307, const double * const W308, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no118_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO118_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W309, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no118_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO118_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W310, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no118_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO118_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W309, const double * const W310, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no119_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO119_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W311, const double * const W312, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no119_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO119_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W312, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no120_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO120_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W313, const double * const W314, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no120_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO120_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W314, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no121_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO121_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W315, const double * const W316, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no121_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO121_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W316, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no122_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO122_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W317, const double * const W318, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no122_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO122_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W318, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no123_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO123_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W319, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no123_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO123_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W320, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no123_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO123_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W319, const double * const W320, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no124_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO124_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W321, const double * const W322, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no124_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO124_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W322, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no125_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO125_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W323, const double * const W324, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no125_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO125_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W324, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no126_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO126_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W325, const double * const W326, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no126_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO126_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W326, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no127_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO127_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W327, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no127_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO127_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W328, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no127_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO127_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W327, const double * const W328, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no128_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO128_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W329, const double * const W330, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no128_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO128_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W330, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no129_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO129_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W331, const double * const W332, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no129_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO129_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W332, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no130_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO130_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W333, const double * const W334, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no130_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO130_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W334, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no131_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO131_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W335, const double * const W336, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no131_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO131_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W336, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no132_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO132_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W337, const double * const W338, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no132_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO132_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W338, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no133_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO133_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W339, const double * const W340, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no133_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO133_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W340, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no134_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO134_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W342, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no134_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO134_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W341, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no134_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO134_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W341, const double * const W342, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no135_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO135_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W343, const double * const W344, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no135_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO135_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W344, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no136_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO136_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W346, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no136_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO136_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W345, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no136_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO136_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W345, const double * const W346, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no137_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO137_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W347, const double * const W348, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no137_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO137_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W348, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no138_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO138_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W349, const double * const W350, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no138_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO138_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W350, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no139_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO139_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W351, const double * const W352, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no139_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO139_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W352, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no140_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO140_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W353, const double * const W354, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no140_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO140_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W354, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no141_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO141_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W355, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no141_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO141_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W355, const double * const W356, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no142_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO142_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W357, const double * const W358, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no142_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO142_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W358, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no143_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO143_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W359, const double * const W360, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no143_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO143_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W360, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no144_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO144_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W361, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no144_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO144_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W361, const double * const W362, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no145_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO145_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W363, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no145_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO145_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W364, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no145_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO145_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W363, const double * const W364, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no146_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO146_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W365, const double * const W366, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no146_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO146_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W366, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no147_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO147_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W367, const double * const W368, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no147_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO147_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W368, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no148_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO148_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W369, const double * const W370, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no148_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO148_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W370, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no149_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO149_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W371, const double * const W372, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no149_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO149_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W372, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no150_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO150_X0_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W373, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no150_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO150_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W374, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no150_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO150_X2_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W373, const double * const W374, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no151_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO151_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W375, const double * const W376, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no151_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO151_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W376, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no152_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO152_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W377, const double * const W378, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no152_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO152_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W378, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no153_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO153_X0_TYPE1_ERI_C)
  (const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W379, const double * const W380, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no153_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO153_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W380, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no154_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO154_X0_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W381, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no154_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO154_X1_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W382, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no154_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO154_X2_TYPE1_ERI_C)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sc0, const FC_INT &ic0, const FC_INT &sj, const FC_INT &ij, 
   const double * const W381, const double * const W382, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no155_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO155_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W383, const double * const W384, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no155_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO155_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W384, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no156_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO156_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W385, const double * const W386, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no156_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO156_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W386, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no157_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO157_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W387, const double * const W388, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no157_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO157_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W388, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no158_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO158_X0_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const V2, const double * const W389, const double * const W390, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no158_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO158_X1_TYPE1_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sw, const FC_INT &iw, 
   const double * const W390, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type2_eri_c,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE2_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W252, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type2_eri_c,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE2_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W274, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type2_eri_c,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE2_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W302, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type2_eri_c,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE2_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W308, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type2_eri_c,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE2_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W356, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type2_eri_c,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE2_ERI_C)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W362, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type0_eri_o,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE0_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W95, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type0_eri_o,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE0_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W96, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type0_eri_o,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE0_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W99, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type0_eri_o,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE0_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const T2, const double * const W100, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type0_eri_o,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE0_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W109, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type0_eri_o,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE0_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W110, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no6_x0_type0_eri_o,G_IF_SIGMA_COOO_COOO_NO6_X0_TYPE0_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const W420, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const V2, const double * const W28, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W29, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no1_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO1_X1_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W29, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const V2, const double * const W30, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const V2, const double * const W31, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W32, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no4_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO4_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W32, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const V2, const double * const W33, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no6_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO6_X0_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W34, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no6_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO6_X1_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W34, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no7_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO7_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const V2, const double * const W35, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no7_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO7_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const W35, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no8_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO8_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W36, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no8_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO8_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W36, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no9_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO9_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const V2, const double * const W39, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no9_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO9_X1_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W39, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no10_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO10_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const V2, const double * const W40, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no10_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO10_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const W40, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no11_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO11_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W41, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no11_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO11_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W41, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no12_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO12_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const V2, const double * const W43, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no12_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO12_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W43, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no13_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO13_X0_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, 
   const double * const T2, const double * const W44, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no13_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO13_X1_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W44, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no14_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO14_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const V2, const double * const W45, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no14_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO14_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const W45, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no15_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO15_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W46, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no15_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO15_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W46, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no16_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO16_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const V2, const double * const W47, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no16_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO16_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const W47, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no17_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO17_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W48, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no17_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO17_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W48, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no18_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO18_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W49, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no18_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO18_X1_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W49, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no19_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO19_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const V2, const double * const W51, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no19_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO19_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W51, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no20_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO20_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const V2, const double * const W52, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no21_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO21_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const V2, const double * const W53, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no22_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO22_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const V2, const double * const W54, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no23_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO23_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W55, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no23_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO23_X1_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W55, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no24_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO24_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const V2, const double * const W56, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no25_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO25_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const V2, const double * const W57, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no26_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO26_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W58, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no26_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO26_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W58, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no27_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO27_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const V2, const double * const W59, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no28_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO28_X0_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W60, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no28_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO28_X1_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W60, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no29_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO29_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const V2, const double * const W61, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no29_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO29_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W61, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no30_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO30_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W62, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no30_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO30_X1_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W62, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no31_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO31_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const V2, const double * const W65, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no31_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO31_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W65, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no32_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO32_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const V2, const double * const W66, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no32_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO32_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W66, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no33_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO33_X0_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W67, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no33_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO33_X1_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W67, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no34_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO34_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const V2, const double * const W69, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no34_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO34_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W69, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no35_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO35_X0_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, 
   const double * const T2, const double * const W70, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no35_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO35_X1_TYPE1_ERI_O)
  (const FC_INT &sa3, const FC_INT &ia3, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W70, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no36_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO36_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const V2, const double * const W71, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no36_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO36_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W71, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no37_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO37_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W72, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no37_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO37_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W72, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no38_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO38_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const V2, const double * const W73, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no38_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO38_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W73, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no39_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO39_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W74, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no39_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO39_X1_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W74, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no40_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO40_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const W75, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no40_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO40_X1_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W75, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no41_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO41_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const V2, const double * const W77, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no41_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO41_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W77, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no42_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO42_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const V2, const double * const W78, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no43_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO43_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const V2, const double * const W79, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no44_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO44_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const V2, const double * const W80, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no45_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO45_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const V2, const double * const W81, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no46_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO46_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const V2, const double * const W82, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no46_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO46_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W82, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no47_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO47_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const V2, const double * const W83, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no47_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO47_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W83, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no48_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO48_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const V2, const double * const W84, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no48_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO48_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W84, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no49_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO49_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const V2, const double * const W85, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no49_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO49_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W85, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no50_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO50_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const V2, const double * const W86, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no51_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO51_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W95, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no52_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO52_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W96, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no53_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO53_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W99, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no54_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO54_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W100, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no55_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO55_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W109, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no56_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO56_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W110, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no57_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO57_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sa4, const FC_INT &ia4, 
   const double * const T2, const double * const W111, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no57_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO57_X1_TYPE1_ERI_O)
  (const FC_INT &sa4, const FC_INT &ia4, 
   const double * const V2, const double * const W111, const double * const W112, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no58_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO58_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sa4, const FC_INT &ia4, 
   const double * const T2, const double * const W125, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no58_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO58_X1_TYPE1_ERI_O)
  (const FC_INT &sa4, const FC_INT &ia4, 
   const double * const V2, const double * const W125, const double * const W126, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no59_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO59_X0_TYPE1_ERI_O)
  (const FC_INT &sa4, const FC_INT &ia4, 
   const double * const V2, const double * const W133, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no59_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO59_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sa4, const FC_INT &ia4, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W134, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no59_x2_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO59_X2_TYPE1_ERI_O)
  (const FC_INT &sa4, const FC_INT &ia4, const FC_INT &sj, const FC_INT &ij, 
   const double * const W133, const double * const W134, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no60_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO60_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sa4, const FC_INT &ia4, 
   const double * const T2, const double * const W145, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no60_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO60_X1_TYPE1_ERI_O)
  (const FC_INT &sa4, const FC_INT &ia4, 
   const double * const V2, const double * const W145, const double * const W146, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no61_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO61_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sa4, const FC_INT &ia4, 
   const double * const T2, const double * const W149, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no61_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO61_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sa4, const FC_INT &ia4, 
   const double * const V2, const double * const W149, const double * const W150, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no62_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO62_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sa4, const FC_INT &ia4, 
   const double * const T2, const double * const W199, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no62_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO62_X1_TYPE1_ERI_O)
  (const FC_INT &sa4, const FC_INT &ia4, 
   const double * const V2, const double * const W199, const double * const W200, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no63_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO63_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sa4, const FC_INT &ia4, 
   const double * const T2, const double * const W203, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no63_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO63_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sa4, const FC_INT &ia4, 
   const double * const V2, const double * const W203, const double * const W204, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no64_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO64_X0_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W391, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no64_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO64_X1_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W391, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no65_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO65_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W392, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no65_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO65_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W392, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no66_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO66_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W393, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no66_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO66_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W393, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no67_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO67_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const V2, const double * const W394, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no67_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO67_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W394, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no68_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO68_X0_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W395, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no68_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO68_X1_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W395, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no69_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO69_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, 
   const double * const V2, const double * const W396, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no69_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO69_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W396, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no70_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO70_X0_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W397, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no70_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO70_X1_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W397, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no71_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO71_X0_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const V2, const double * const W398, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no71_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO71_X1_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W398, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no72_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO72_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const V2, const double * const W399, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no72_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO72_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W399, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no73_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO73_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const V2, const double * const W400, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no73_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO73_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W400, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no74_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO74_X0_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const V2, const double * const W401, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no74_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO74_X1_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W401, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no75_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO75_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W402, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no75_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO75_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W402, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no76_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO76_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const V2, const double * const W416, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no76_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO76_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W416, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no77_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO77_X0_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, 
   const double * const V2, const double * const W417, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no77_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO77_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W417, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no78_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO78_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W418, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no78_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO78_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W418, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no79_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO79_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W419, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no79_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO79_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W419, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no80_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO80_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W420, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no81_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO81_X0_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, const FC_INT &sk, const FC_INT &ik, 
   const double * const V2, const double * const W421, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no81_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO81_X1_TYPE1_ERI_O)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, const FC_INT &sk, const FC_INT &ik, 
   const double * const T2, const double * const W421, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no82_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO82_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W422, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no82_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO82_X1_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W422, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no83_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO83_X0_TYPE1_ERI_O)
  (const FC_INT &sk, const FC_INT &ik, 
   const double * const V2, const double * const W423, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no83_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO83_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sk, const FC_INT &ik, 
   const double * const T2, const double * const W423, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no84_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO84_X0_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sk, const FC_INT &ik, 
   const double * const V2, const double * const W424, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no84_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO84_X1_TYPE1_ERI_O)
  (const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, const FC_INT &sk, const FC_INT &ik, 
   const double * const T2, const double * const W424, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no85_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO85_X0_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W425, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no85_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO85_X1_TYPE1_ERI_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W425, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no86_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO86_X0_TYPE1_ERI_O)
  (const FC_INT &sk, const FC_INT &ik, 
   const double * const V2, const double * const W426, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no86_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO86_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sk, const FC_INT &ik, 
   const double * const T2, const double * const W426, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no87_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO87_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W427, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no87_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO87_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W427, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no88_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO88_X0_TYPE1_ERI_O)
  (const FC_INT &sk, const FC_INT &ik, 
   const double * const V2, const double * const W428, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no88_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO88_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sk, const FC_INT &ik, 
   const double * const T2, const double * const W428, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no89_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO89_X0_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const V2, const double * const W429, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no89_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO89_X1_TYPE1_ERI_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const T2, const double * const W429, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no90_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO90_X0_TYPE1_ERI_O)
  (const FC_INT &sk, const FC_INT &ik, 
   const double * const V2, const double * const W430, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no90_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO90_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sk, const FC_INT &ik, 
   const double * const T2, const double * const W430, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no91_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO91_X0_TYPE1_ERI_O)
  (const FC_INT &sk, const FC_INT &ik, 
   const double * const V2, const double * const W431, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no91_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO91_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sk, const FC_INT &ik, 
   const double * const T2, const double * const W431, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no92_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO92_X0_TYPE1_ERI_O)
  (const FC_INT &sk, const FC_INT &ik, 
   const double * const V2, const double * const W432, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no92_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO92_X1_TYPE1_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, const FC_INT &sk, const FC_INT &ik, 
   const double * const T2, const double * const W432, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W28, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W30, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W31, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W33, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W52, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W53, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no6_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO6_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W54, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no7_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO7_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W56, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no8_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO8_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W57, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no9_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO9_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W59, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no10_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO10_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W78, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no11_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO11_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W79, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no12_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO12_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W80, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no13_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO13_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W81, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no14_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO14_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W86, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no15_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO15_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W112, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no16_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO16_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W126, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no17_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO17_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W146, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no18_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO18_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W150, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no19_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO19_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W200, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no20_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO20_X0_TYPE2_ERI_O)
  (const FC_INT &sj, const FC_INT &ij, 
   const double * const W204, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type1_d4c_o,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE1_D4C_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sa3, const FC_INT &ia3, const FC_INT &sj, const FC_INT &ij, 
   const double * const C5, const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type1_d4c_o,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE1_D4C_O)
  (const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const C5, const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type1_d4c_o,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE1_D4C_O)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sa2, const FC_INT &ia2, const FC_INT &sj, const FC_INT &ij, 
   const double * const C5, const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type1_d4c_o,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE1_D4C_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &si, const FC_INT &ii, const FC_INT &sj, const FC_INT &ij, 
   const double * const C5, const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type1_d4c_o,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE1_D4C_O)
  (const FC_INT &sa1, const FC_INT &ia1, const FC_INT &sj, const FC_INT &ij, 
   const double * const C5, const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type1_d4c_o,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE1_D4C_O)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sj, const FC_INT &ij, const FC_INT &sk, const FC_INT &ik, 
   const double * const C5, const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 