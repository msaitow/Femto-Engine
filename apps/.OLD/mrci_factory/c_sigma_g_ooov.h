

#ifndef C_SIGMA_G_OOOV_H
#define C_SIGMA_G_OOOV_H

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

void FC_FUNC(g_if_sigma_g_ooov_no0_x0, G_IF_SIGMA_G_OOOV_NO0_X0)
  (const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_g_ooov_no1_x0, G_IF_SIGMA_G_OOOV_NO1_X0)
  (const double * const X, const double * const S0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_g_ooov_y1, G_IF_SIGMA_G_OOOV_Y1)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_g_ooov_no0_x1, G_IF_SIGMA_G_OOOV_NO0_X1)
  (const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const Y0, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_g_ooov_no1_x1, G_IF_SIGMA_G_OOOV_NO1_X1)
  (const double * const X, const double * const S0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_g_ooov_y2, G_IF_SIGMA_G_OOOV_Y2)
  (const FC_INT &sc1, const FC_INT &ic1, 
   const double * const V2, const double * const Y1, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_g_ooov_no0_x2, G_IF_SIGMA_G_OOOV_NO0_X2)
  (const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const Y1, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_g_ooov_no1_x2, G_IF_SIGMA_G_OOOV_NO1_X2)
  (const double * const X, const double * const S0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_g_ooov_no0_x3, G_IF_SIGMA_G_OOOV_NO0_X3)
  (const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_g_ooov_no1_x3, G_IF_SIGMA_G_OOOV_NO1_X3)
  (const double * const X, const double * const S0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_g_ooov_no0_x4, G_IF_SIGMA_G_OOOV_NO0_X4)
  (const FC_INT &sv1, const FC_INT &iv1, 
   const double * const T2, const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_g_ooov_no1_x4, G_IF_SIGMA_G_OOOV_NO1_X4)
  (const double * const X, const double * const S0, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

      
 }     
       
       
 #endif
       
       
 