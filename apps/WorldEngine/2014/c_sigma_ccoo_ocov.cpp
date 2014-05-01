                                                                                
#include <orz/orz.h>
#include <orz/openmp.h>
#include <orz/cblas.h>
#include <orz/clapack.h>
#include <tensor/tensor.h>
#include <sci/hint/para_disttools.h>
#include <sci/icmr/mr.h>
#include <sci/icmr/mr_f.h>
#include <sci/icmr/mrclass_input.h>
#include <sci/icmr/mrclass_symblock.h>
#include <sci/icmr/mrclass_hintmo.h>
#include <sci/icmr/mrclass_rdmpack.h>
#include <sci/icmr/mrclass_bareamppack.h>
#include <sci/icmr/mrclass_orthamppack.h>
#include <sci/icmr/diaghessian.h>
#include <sci/icmr/symamp2.h>
#include <sci/icmr/femto/femto.h>
#include <sci/icmr/Femto/elems/c_sigma_ccoo_ocov.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//  8 8888888888   8 8888888888            ,8.       ,8.    8888888 8888888888 ,o888888o.     
//  8 8888         8 8888                 ,888.     ,888.         8 8888    . 8888     `88.   
//  8 8888         8 8888                .`8888.   .`8888.        8 8888   ,8 8888       `8b  
//  8 8888         8 8888               ,8.`8888. ,8.`8888.       8 8888   88 8888        `8b 
//  8 888888888888 8 888888888888      ,8'8.`8888,8^8.`8888.      8 8888   88 8888         88 
//  8 8888         8 8888             ,8' `8.`8888' `8.`8888.     8 8888   88 8888         88 
//  8 8888         8 8888            ,8'   `8.`88'   `8.`8888.    8 8888   88 8888        ,8P 
//  8 8888         8 8888           ,8'     `8.`'     `8.`8888.   8 8888   `8 8888       ,8P  
//  8 8888         8 8888          ,8'       `8        `8.`8888.  8 8888    ` 8888     ,88'   
//  8 8888         8 888888888888 ,8'         `         `8.`8888. 8 8888       `8888888P'     

//                                   Generated date : Sun Apr 20 10:26:20 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccoo_ocov(const orz::mr::Input &ctinp,                                    
                                  const orz::mr::SymBlockInfo &symblockinfo,                                 
                                  const orz::mr::HintMO &hintmo,                                             
                                  const int alloc_type,                                                      
                                  const orz::mr::BareAmpPack &T2,                             
                                  const int num_sigma) {
                                                                                                                 
                                                                                                                 
  // set up nmo nclosed, nocc                                                                                    
  const FC_INT nclosed = ctinp.nclosed();                                                                        
  const FC_INT nocc    = ctinp.nocc();                                                                           
  const FC_INT nvir    = ctinp.nvir();                                                                           
  const FC_INT nmo     = nclosed + nocc + nvir;                                                                  
  const FC_INT nir     = symblockinfo.nir();                                                                     
  const FC_INT * const nsym    = symblockinfo.nsym().cptr();                                                     
  const FC_INT * const psym    = symblockinfo.psym().cptr();                                                     
  const FC_INT * const amo2imo = symblockinfo.amo2imo().cptr();                                                  
                                                                                                                 
  std::ostringstream stm;                                                                                        
  stm << num_sigma;                                                                                              
  std::string name_of_sigma = "S2" + stm.str() + "]"; // Name of the Sigma vector  
  orz::mr::BareAmpPack retval                                                                                    
    = orz::mr::BareAmpPack(ctinp, symblockinfo, name_of_sigma, alloc_type); // Sigma(a, a', e, e') tensor        
                                                                                                                 
  orz::DTensor S2b; // Container of S2_aae,[b] tensor                                   
                                                                                                                 
  orz::DTensor T2b; // Container of T2_aae,[b] tensor                                             
  // set nproc, myrank                      
  const int nproc = orz::world().size();    
  const int myrank = orz::world().rank();   

  orz::DTensor moint1 = hintmo.int1(); // Setting up one-body integrals                                         
  const orz::DTensor moint1_sym = (myrank == 0) ? orz::mr::sympack_int1(symblockinfo, moint1) : orz::DTensor(); // moint1=(IR-COV index)
  orz::DTensor V2(nmo,nmo,nmo);                                                                    
  double * const V2_ptr = V2.cptr();                                                  

  // Timing object
  orz::ProgressTimer time_sigma(false);

//-@ERI.contractions(begin)

//-@loadERI(c,begin)
  //*-- FEMTO begins --//*
  // Label : eri_c
  {


  for(int s_eri = 0;s_eri < nir;++s_eri){ 
  for(int i_eri = symblockinfo.psym()(s_eri,I_C,I_BEGIN);i_eri <= symblockinfo.psym()(s_eri,I_C,I_END);++i_eri){ 
  if(hintmo.iproc_havingimo()[i_eri] == myrank) {           
  // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
  V2 <<= 0.0;                                                                          
  shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(i_eri);
  for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
    // Load a signle record of integals                                                             
    const int &imo2 = loadbuf_ptr->i0;                                                              
    const int &imo3 = loadbuf_ptr->i1;                                                              
    const int &imo4 = loadbuf_ptr->i2;                                                              
    const double &v = loadbuf_ptr->v;                                                               
    V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
  }                                                                                                 
  const orz::DTensor V2_sym = orz::mr::sympack_int2(symblockinfo, i_eri, s_eri, V2); // V2=(IR-COV index) 

  //*-- Entering to take the type 1 contractions --*//
//-@type(1).contraction(begin)
  // -- Title : sigma_ccoo_ocov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W0(w,x,a1,a0,a3,a2) += (    1.00000000) V2(x,a3,v0,a2) T2(w,a1,v0,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) D3(j,a3,i,a2,a1,a0) W0(w,x,a1,a0,a3,a2) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W0caaa_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, sx^sa0));
    FC_FUNC(g_if_sigma_ccoo_ocov_no0_x0_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO0_X0_TYPE1_ERI_C)
      (sa0, ia0, sx, ix, T2b.cptr(), V2_sym.cptr(), W0caaa_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ocov_no0_x1_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO0_X1_TYPE1_ERI_C)
        (sa0, ia0, sj, ij, sx, ix, W0caaa_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(x,w,a1,a0,a3,a2) += (    1.00000000) V2(w,a3,v0,a2) T2(x,a1,v0,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) D3(j,a2,i,a3,a1,a0) W1(x,w,a1,a0,a3,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W1caaa_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, sw^sa0));
    FC_FUNC(g_if_sigma_ccoo_ocov_no1_x0_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO1_X0_TYPE1_ERI_C)
      (sa0, ia0, sw, iw, T2b.cptr(), V2_sym.cptr(), W1caaa_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ocov_no1_x1_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO1_X1_TYPE1_ERI_C)
        (sa0, ia0, sj, ij, sw, iw, W1caaa_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W10(x,i,a0,v0) += (    1.00000000) V2(x,a2,v0,a1) D2(i,a1,a0,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,a0,v0,j) W10(x,i,a0,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W10aav_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ocov_no2_x0_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO2_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W10aav_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ocov_no2_x1_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO2_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W10aav_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W11(w,i,a0,v0) += (    1.00000000) V2(w,a2,v0,a1) D2(i,a2,a0,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a0,v0,j) W11(w,i,a0,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W11aav_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ocov_no3_x0_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO3_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W11aav_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ocov_no3_x1_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO3_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W11aav_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W14(x,a0,i,v0) += (    1.00000000) V2(x,i,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,a0,v0,j) W14(x,a0,i,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W14aav_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ocov_no4_x0_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO4_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W14aav_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ocov_no4_x1_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO4_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W14aav_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W15(w,a0,i,v0) += (    1.00000000) V2(w,i,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,a0,v0,j) W15(w,a0,i,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W15aav_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ocov_no5_x0_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO5_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W15aav_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ocov_no5_x1_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO5_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W15aav_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W24(w,a0,i,v0) += (    1.00000000) V2(w,a1,v0,i) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a0,v0,j) W24(w,a0,i,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W24aav_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ocov_no6_x0_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO6_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W24aav_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ocov_no6_x1_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO6_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W24aav_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W25(x,a0,i,v0) += (    1.00000000) V2(x,a1,v0,i) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,a0,v0,j) W25(x,a0,i,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W25aav_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ocov_no7_x0_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO7_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W25aav_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ocov_no7_x1_type1_eri_c,G_IF_SIGMA_CCOO_OCOV_NO7_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W25aav_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

//-@type(1).contraction(end)

  } // End myrank
  }
  }
  orz::world().barrier();

  } // End Femto
  //*-- FEMTO ends --//*

//-@loadERI(c,end)

//-@loadERI(a,begin)
  //*-- FEMTO begins --//*
  // Label : eri_o
  {


  for(int s_eri = 0;s_eri < nir;++s_eri){ 
  for(int i_eri = symblockinfo.psym()(s_eri,I_O,I_BEGIN);i_eri <= symblockinfo.psym()(s_eri,I_O,I_END);++i_eri){ 
  if(hintmo.iproc_havingimo()[i_eri] == myrank) {           
  // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
  V2 <<= 0.0;                                                                          
  shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(i_eri);
  for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
    // Load a signle record of integals                                                             
    const int &imo2 = loadbuf_ptr->i0;                                                              
    const int &imo3 = loadbuf_ptr->i1;                                                              
    const int &imo4 = loadbuf_ptr->i2;                                                              
    const double &v = loadbuf_ptr->v;                                                               
    V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
  }                                                                                                 
  const orz::DTensor V2_sym = orz::mr::sympack_int2(symblockinfo, i_eri, s_eri, V2); // V2=(IR-COV index) 

  //*-- Entering to take the type 1 contractions --*//
//-@type(1).contraction(begin)
  // -- Title : sigma_ccoo_ocov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W18(x,a0,j,v0) += (    1.00000000) V2(j,x,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,a0,v0,i) W18(x,a0,j,v0) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W18cav_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ccoo_ocov_no0_x0_type1_eri_o,G_IF_SIGMA_CCOO_OCOV_NO0_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W18cav_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_ocov_no0_x1_type1_eri_o,G_IF_SIGMA_CCOO_OCOV_NO0_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W18cav_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W19(w,a0,j,v0) += (    1.00000000) V2(j,w,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a0,v0,i) W19(w,a0,j,v0) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W19cav_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ccoo_ocov_no1_x0_type1_eri_o,G_IF_SIGMA_CCOO_OCOV_NO1_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W19cav_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_ocov_no1_x1_type1_eri_o,G_IF_SIGMA_CCOO_OCOV_NO1_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W19cav_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

//-@type(1).contraction(end)

  } // End myrank
  }
  }
  orz::world().barrier();

  } // End Femto
  //*-- FEMTO ends --//*

//-@loadERI(a,end)

//-@loadERI(v,begin)
  //*-- FEMTO begins --//*
  // Label : eri_v
  {


  for(int s_eri = 0;s_eri < nir;++s_eri){ 
  for(int i_eri = symblockinfo.psym()(s_eri,I_V,I_BEGIN);i_eri <= symblockinfo.psym()(s_eri,I_V,I_END);++i_eri){ 
  if(hintmo.iproc_havingimo()[i_eri] == myrank) {           
  // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
  V2 <<= 0.0;                                                                          
  shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(i_eri);
  for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
    // Load a signle record of integals                                                             
    const int &imo2 = loadbuf_ptr->i0;                                                              
    const int &imo3 = loadbuf_ptr->i1;                                                              
    const int &imo4 = loadbuf_ptr->i2;                                                              
    const double &v = loadbuf_ptr->v;                                                               
    V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
  }                                                                                                 
  const orz::DTensor V2_sym = orz::mr::sympack_int2(symblockinfo, i_eri, s_eri, V2); // V2=(IR-COV index) 

  //*-- Entering to take the type 1 contractions --*//
//-@type(1).contraction(begin)
  // -- Title : sigma_ccoo_ocov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W2(w,j,a2,v0) += (    1.00000000) T2(a1,w,a0,v0) D2(j,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(v0,a2,x,i) W2(w,j,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W2ca_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ocov_no0_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO0_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W2ca_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ocov_no0_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO0_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W2ca_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W3(x,j,a2,v0) += (    1.00000000) T2(a1,x,a0,v0) D2(j,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(v0,a2,w,i) W3(x,j,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W3ca_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ocov_no1_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO1_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W3ca_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ocov_no1_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO1_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W3ca_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W4(x,j,a0,v0) += (    1.00000000) V2(v0,a1,x,a2) D2(j,a2,a0,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(a0,w,i,v0) W4(x,j,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W4ca_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ocov_no2_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO2_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W4ca_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ocov_no2_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO2_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W4ca_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W5(w,j,a0,v0) += (    1.00000000) V2(v0,a1,w,a2) D2(j,a1,a0,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(a0,x,i,v0) W5(w,j,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W5ca_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ocov_no3_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO3_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W5ca_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ocov_no3_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO3_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W5ca_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W6(x,j,a2,v0) += (    1.00000000) T2(a1,x,a0,v0) D2(j,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(v0,i,w,a2) W6(x,j,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W6ca_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ocov_no4_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO4_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W6ca_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ocov_no4_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO4_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W6ca_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W7(w,j,a2,v0) += (    1.00000000) T2(a1,w,a0,v0) D2(j,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(v0,i,x,a2) W7(w,j,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W7ca_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ocov_no5_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO5_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W7ca_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ocov_no5_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO5_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W7ca_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W8(w,i,a2,v0) += (    1.00000000) T2(a1,w,a0,v0) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(v0,a2,x,j) W8(w,i,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W8caa_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccoo_ocov_no6_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO6_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W8caa_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ocov_no6_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO6_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W8caa_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W9(x,i,a2,v0) += (    1.00000000) T2(a1,x,a0,v0) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(v0,a2,w,j) W9(x,i,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W9caa_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccoo_ocov_no7_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO7_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W9caa_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ocov_no7_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO7_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W9caa_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W12(x,i,a2,v0) += (    1.00000000) T2(a1,x,a0,v0) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(v0,j,w,a2) W12(x,i,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W12caa_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccoo_ocov_no8_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO8_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W12caa_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ocov_no8_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO8_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W12caa_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W13(w,i,a2,v0) += (    1.00000000) T2(a1,w,a0,v0) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(v0,j,x,a2) W13(w,i,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W13caa_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccoo_ocov_no9_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO9_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W13caa_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ocov_no9_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO9_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W13caa_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W16(x,v0) += (    1.00000000) T2(a1,x,a0,v0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) V2(v0,j,w,i) W16(x,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W16c_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xc(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccoo_ocov_no10_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO10_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W16c_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ocov_no10_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO10_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W16c_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W17(w,v0) += (    1.00000000) T2(a1,w,a0,v0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(v0,j,x,i) W17(w,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W17c_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xc(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccoo_ocov_no11_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO11_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W17c_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ocov_no11_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO11_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W17c_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W20(w,a0,j,v0) += (    1.00000000) V2(v0,j,w,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(a0,x,i,v0) W20(w,a0,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W20ca_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ocov_no12_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO12_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W20ca_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ocov_no12_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO12_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W20ca_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W21(x,a0,j,v0) += (    1.00000000) V2(v0,j,x,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(a0,w,i,v0) W21(x,a0,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W21ca_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ocov_no13_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO13_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W21ca_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ocov_no13_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO13_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W21ca_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W22(w,v0) += (    1.00000000) T2(a1,w,a0,v0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) V2(v0,i,x,j) W22(w,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W22c_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xc(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccoo_ocov_no14_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO14_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W22c_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ocov_no14_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO14_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W22c_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W23(x,v0) += (    1.00000000) T2(a1,x,a0,v0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(v0,i,w,j) W23(x,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W23c_sigma_ccoo_ocov(orz::mr::sizeof_sympack_Xc(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccoo_ocov_no15_x0_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO15_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W23c_sigma_ccoo_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ocov_no15_x1_type1_eri_v,G_IF_SIGMA_CCOO_OCOV_NO15_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W23c_sigma_ccoo_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

//-@type(1).contraction(end)

  } // End myrank
  }
  }
  orz::world().barrier();

  } // End Femto
  //*-- FEMTO ends --//*

//-@loadERI(v,end)

//-@ERI.contractions(end)

//-@D4C.contractions(begin)

//-@D4C.contractions(end)

  // Do timing!
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccoo_ocov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
