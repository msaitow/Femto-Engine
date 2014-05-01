

#ifndef C_SIGMA_OOOV_COOO_H
#define C_SIGMA_OOOV_COOO_H

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

void FC_FUNC(g_if_sigma_ooov_cooo_no0_x0_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO0_X0_TYPE1_ERI_C)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W12, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no0_x1_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO0_X1_TYPE1_ERI_C)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const W12, const double * const W13, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no0_x2_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO0_X2_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W13, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no1_x0_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO1_X0_TYPE1_ERI_C)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W20, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no1_x1_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO1_X1_TYPE1_ERI_C)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const W20, const double * const W21, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no1_x2_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO1_X2_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W21, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no2_x0_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO2_X0_TYPE1_ERI_C)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W32, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no2_x1_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO2_X1_TYPE1_ERI_C)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const W32, const double * const W33, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no2_x2_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO2_X2_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W33, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no3_x0_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO3_X0_TYPE1_ERI_C)
  (const FC_INT &sa0, const FC_INT &ia0, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const T2, const double * const W40, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no3_x1_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO3_X1_TYPE1_ERI_C)
  (const FC_INT &sc0, const FC_INT &ic0, 
   const double * const W40, const double * const W41, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no3_x2_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO3_X2_TYPE1_ERI_C)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sc0, const FC_INT &ic0, 
   const double * const V2, const double * const W41, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no0_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO0_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W10, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no1_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO1_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W15, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no2_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO2_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W17, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no3_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO3_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W19, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no4_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO4_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W23, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no5_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO5_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W24, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no6_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO6_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W28, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no7_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO7_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W34, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no8_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO8_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W42, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no9_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO9_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W44, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no10_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO10_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W46, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no11_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO11_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W48, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no12_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO12_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W50, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no13_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO13_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W52, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no14_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO14_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W54, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no15_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO15_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W56, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no16_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO16_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W58, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no17_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO17_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W60, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no18_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO18_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W62, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no19_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO19_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W65, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no20_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO20_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W67, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no21_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO21_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W68, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no22_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO22_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W70, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no23_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO23_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W72, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no24_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO24_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W74, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no25_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO25_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W76, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no26_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO26_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W79, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no27_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO27_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W80, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no28_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO28_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W82, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no29_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO29_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W84, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no30_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO30_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W86, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no31_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO31_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W89, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no32_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO32_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W91, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no33_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO33_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W93, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no34_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO34_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W95, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no35_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO35_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W96, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no36_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO36_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W98, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no37_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO37_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W100, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no38_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO38_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W102, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no39_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO39_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W104, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no40_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO40_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W106, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no41_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO41_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W108, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no42_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO42_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W110, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no43_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO43_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W112, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no44_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO44_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W114, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no45_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO45_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W116, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no46_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO46_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W119, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no47_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO47_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W121, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no48_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO48_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W122, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no49_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO49_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W124, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no50_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO50_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W126, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no51_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO51_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W128, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no52_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO52_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W130, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no53_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO53_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W133, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no54_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO54_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W134, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no55_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO55_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W136, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no56_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO56_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W138, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no57_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO57_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W140, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no58_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO58_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W143, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no59_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO59_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W145, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no60_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO60_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W147, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no61_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO61_X0_TYPE0_ERI_V)
  (const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W149, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no62_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO62_X0_TYPE0_ERI_V)
  (const FC_INT &sk, const FC_INT &ik, 
   const double * const T2, const double * const W150, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no63_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO63_X0_TYPE0_ERI_V)
  (const FC_INT &sk, const FC_INT &ik, 
   const double * const T2, const double * const W151, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no0_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO0_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const V2, const double * const W0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no0_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO0_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const W0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no1_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO1_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const V2, const double * const W1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no1_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO1_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const W1, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no2_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO2_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sk, const FC_INT &ik, 
   const double * const T2, const double * const V2, const double * const W2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no2_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO2_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sk, const FC_INT &ik, 
   const double * const W2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no3_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO3_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sk, const FC_INT &ik, 
   const double * const T2, const double * const V2, const double * const W3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no3_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO3_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sk, const FC_INT &ik, 
   const double * const W3, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no4_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO4_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const V2, const double * const W4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no4_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO4_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W4, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no5_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO5_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa3, const FC_INT &ia3, 
   const double * const T2, const double * const V2, const double * const W5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no5_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO5_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W5, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no6_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO6_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const V2, const double * const W6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no6_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO6_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W6, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no7_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO7_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const V2, const double * const W7, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no7_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO7_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W7, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no8_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO8_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa1, const FC_INT &ia1, 
   const double * const T2, const double * const V2, const double * const W8, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no8_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO8_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W8, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no9_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO9_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa2, const FC_INT &ia2, 
   const double * const T2, const double * const V2, const double * const W9, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no9_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO9_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W9, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no10_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO10_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W10, const double * const W11, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no10_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO10_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W11, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no11_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO11_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W14, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no11_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO11_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W14, const double * const W15, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no12_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO12_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W16, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no12_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO12_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W16, const double * const W17, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no13_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO13_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W18, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no13_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO13_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W18, const double * const W19, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no14_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO14_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W22, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no14_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO14_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W22, const double * const W23, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no15_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO15_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W24, const double * const W25, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no15_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO15_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W25, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no16_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO16_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W26, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no16_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO16_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W26, const double * const W27, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no16_x2_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO16_X2_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const W27, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no17_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO17_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W28, const double * const W29, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no17_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO17_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W29, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no18_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO18_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W30, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no18_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO18_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W30, const double * const W31, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no18_x2_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO18_X2_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const W31, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no19_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO19_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W34, const double * const W35, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no19_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO19_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W35, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no20_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO20_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const V2, const double * const W36, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no20_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO20_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W36, const double * const W37, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no20_x2_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO20_X2_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W37, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no21_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO21_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const V2, const double * const W38, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no21_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO21_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sa0, const FC_INT &ia0, 
   const double * const T2, const double * const W38, const double * const W39, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no21_x2_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO21_X2_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W39, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no22_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO22_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W42, const double * const W43, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no22_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO22_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W43, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no23_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO23_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W44, const double * const W45, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no23_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO23_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W45, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no24_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO24_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W46, const double * const W47, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no24_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO24_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W47, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no25_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO25_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W48, const double * const W49, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no25_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO25_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W49, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no26_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO26_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W50, const double * const W51, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no26_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO26_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W51, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no27_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO27_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W52, const double * const W53, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no27_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO27_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W53, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no28_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO28_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W54, const double * const W55, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no28_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO28_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W55, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no29_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO29_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W56, const double * const W57, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no29_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO29_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W57, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no30_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO30_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W58, const double * const W59, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no30_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO30_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W59, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no31_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO31_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W60, const double * const W61, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no31_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO31_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W61, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no32_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO32_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W62, const double * const W63, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no32_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO32_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W63, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no33_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO33_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W64, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no33_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO33_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W64, const double * const W65, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no34_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO34_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W66, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no34_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO34_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W66, const double * const W67, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no35_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO35_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W68, const double * const W69, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no35_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO35_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W69, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no36_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO36_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W70, const double * const W71, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no36_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO36_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W71, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no37_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO37_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W72, const double * const W73, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no37_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO37_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W73, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no38_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO38_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W74, const double * const W75, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no38_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO38_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W75, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no39_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO39_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W76, const double * const W77, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no39_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO39_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W77, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no40_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO40_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W78, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no40_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO40_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W78, const double * const W79, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no41_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO41_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W80, const double * const W81, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no41_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO41_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W81, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no42_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO42_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W82, const double * const W83, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no42_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO42_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W83, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no43_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO43_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W84, const double * const W85, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no43_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO43_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W85, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no44_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO44_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W86, const double * const W87, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no44_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO44_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W87, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no45_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO45_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W88, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no45_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO45_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W88, const double * const W89, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no46_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO46_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W90, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no46_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO46_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W90, const double * const W91, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no47_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO47_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W92, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no47_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO47_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W92, const double * const W93, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no48_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO48_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W94, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no48_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO48_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W94, const double * const W95, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no49_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO49_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W96, const double * const W97, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no49_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO49_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W97, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no50_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO50_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W98, const double * const W99, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no50_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO50_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W99, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no51_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO51_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W100, const double * const W101, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no51_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO51_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W101, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no52_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO52_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W102, const double * const W103, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no52_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO52_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W103, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no53_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO53_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W104, const double * const W105, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no53_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO53_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W105, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no54_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO54_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W106, const double * const W107, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no54_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO54_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W107, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no55_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO55_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W108, const double * const W109, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no55_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO55_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W109, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no56_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO56_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W110, const double * const W111, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no56_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO56_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W111, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no57_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO57_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W112, const double * const W113, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no57_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO57_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W113, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no58_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO58_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W114, const double * const W115, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no58_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO58_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W115, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no59_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO59_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W116, const double * const W117, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no59_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO59_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W117, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no60_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO60_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W118, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no60_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO60_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W118, const double * const W119, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no61_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO61_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W120, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no61_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO61_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W120, const double * const W121, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no62_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO62_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W122, const double * const W123, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no62_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO62_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W123, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no63_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO63_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W124, const double * const W125, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no63_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO63_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W125, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no64_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO64_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W126, const double * const W127, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no64_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO64_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W127, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no65_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO65_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W128, const double * const W129, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no65_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO65_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W129, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no66_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO66_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W130, const double * const W131, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no66_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO66_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W131, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no67_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO67_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W132, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no67_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO67_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W132, const double * const W133, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no68_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO68_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W134, const double * const W135, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no68_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO68_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W135, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no69_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO69_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W136, const double * const W137, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no69_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO69_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W137, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no70_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO70_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W138, const double * const W139, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no70_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO70_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W139, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no71_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO71_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W140, const double * const W141, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no71_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO71_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W141, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no72_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO72_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W142, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no72_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO72_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W142, const double * const W143, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no73_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO73_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W144, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no73_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO73_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W144, const double * const W145, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no74_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO74_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W146, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no74_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO74_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W146, const double * const W147, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no75_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO75_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W148, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no75_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO75_X1_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W148, const double * const W149, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no76_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO76_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W150, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_ooov_cooo_no77_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO77_X0_TYPE1_ERI_V)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const W151, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 