

#ifndef C_SIGMA_OOVV_OOVV_H
#define C_SIGMA_OOVV_OOVV_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//    o__ __o__/_                            o                         
//   <|    v                                <|>                        
//   < >                                    < >                        
//    |         o__  __o   \o__ __o__ __o    |        o__ __o         
//    o__/_    /v      |>   |     |     |>   o__/_   /v     v\        
//    |       />      //   / \   / \   / \   |      />       <\    
//   <o>      \o    o/     \o/   \o/   \o/   |      \         /   
//    |        v\  /v __o   |     |     |    o       o       o        
//   / \        <\/> __/>  / \   / \   / \   <\__    <\__ __/>  

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x0_type0_noeri,G_IF_SIGMA_OOVV_OOVV_NO0_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const Fc0, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x0_type0_noeri,G_IF_SIGMA_OOVV_OOVV_NO1_X0_TYPE0_NOERI)
  (const double * const W0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x1_type0_noeri,G_IF_SIGMA_OOVV_OOVV_NO1_X1_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W0, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no2_x0_type0_noeri,G_IF_SIGMA_OOVV_OOVV_NO2_X0_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const W1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no2_x1_type0_noeri,G_IF_SIGMA_OOVV_OOVV_NO2_X1_TYPE0_NOERI)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const W1, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no3_x0_type0_noeri,G_IF_SIGMA_OOVV_OOVV_NO3_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const W2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no3_x1_type0_noeri,G_IF_SIGMA_OOVV_OOVV_NO3_X1_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no4_x0_type0_noeri,G_IF_SIGMA_OOVV_OOVV_NO4_X0_TYPE0_NOERI)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x0_type1_eri_v,G_IF_SIGMA_OOVV_OOVV_NO0_X0_TYPE1_ERI_V)
  (const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W3, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x1_type1_eri_v,G_IF_SIGMA_OOVV_OOVV_NO0_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W3, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x0_type1_eri_v,G_IF_SIGMA_OOVV_OOVV_NO1_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W4, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x1_type1_eri_v,G_IF_SIGMA_OOVV_OOVV_NO1_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W4, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no2_x0_type1_eri_v,G_IF_SIGMA_OOVV_OOVV_NO2_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W5, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no2_x1_type1_eri_v,G_IF_SIGMA_OOVV_OOVV_NO2_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W5, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no3_x0_type1_eri_v,G_IF_SIGMA_OOVV_OOVV_NO3_X0_TYPE1_ERI_V)
  (const FC_INT &sv0, const FC_INT &iv0, 
   const double * const V2, const double * const W6, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no3_x1_type1_eri_v,G_IF_SIGMA_OOVV_OOVV_NO3_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const W6, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no4_x0_type1_eri_v,G_IF_SIGMA_OOVV_OOVV_NO4_X0_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sv0, const FC_INT &iv0, 
   const double * const T2, const double * const V2, const double * const W7, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

void FC_FUNC(g_if_sigma_oovv_oovv_no4_x1_type1_eri_v,G_IF_SIGMA_OOVV_OOVV_NO4_X1_TYPE1_ERI_V)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const W7, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);

      
 }     
       
       
 #endif
       
       
 