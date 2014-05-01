                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ccov_coov.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//      ______                  __           
//     / ____/___   ____ ___   / /_ ____     
//    / /_   / _ \ / __ `__ \ / __// __ \ 
//   / __/  /  __// / / / / // /_ / /_/ /    
//  /_/     \___//_/ /_/ /_/ \__/ \____/  

//                                   Generated date : Sun Apr 20 10:26:23 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccov_coov(const orz::mr::Input &ctinp,                                    
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

  //*-- FEMTO begins --//*
  // Label : noeri
  {


  //*-- Entering to take the type 0 contractions --*//

  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W0(x,i,a2,a) += (    1.00000000) T2(x,a1,a0,a) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) Fc1(w,a2) W0(x,i,a2,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W0caa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no0_x0_type0_noeri,G_IF_SIGMA_CCOV_COOV_NO0_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W0caa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_coov_no0_x1_type0_noeri,G_IF_SIGMA_CCOV_COOV_NO0_X1_TYPE0_NOERI)
      (sa, ia, W0caa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(w,i,a2,a) += (    1.00000000) T2(w,a1,a0,a) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) Fc1(x,a2) W1(w,i,a2,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W1caa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no1_x0_type0_noeri,G_IF_SIGMA_CCOV_COOV_NO1_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W1caa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_coov_no1_x1_type0_noeri,G_IF_SIGMA_CCOV_COOV_NO1_X1_TYPE0_NOERI)
      (sa, ia, W1caa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(x,a) += (    1.00000000) T2(x,a1,a0,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) Fc1(w,i) W2(x,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W2c_sigma_ccov_coov(orz::mr::sizeof_sympack_Xc(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no2_x0_type0_noeri,G_IF_SIGMA_CCOV_COOV_NO2_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W2c_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_coov_no2_x1_type0_noeri,G_IF_SIGMA_CCOV_COOV_NO2_X1_TYPE0_NOERI)
      (sa, ia, W2c_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(w,a) += (    1.00000000) T2(w,a1,a0,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) Fc1(x,i) W3(w,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W3c_sigma_ccov_coov(orz::mr::sizeof_sympack_Xc(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no3_x0_type0_noeri,G_IF_SIGMA_CCOV_COOV_NO3_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W3c_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_coov_no3_x1_type0_noeri,G_IF_SIGMA_CCOV_COOV_NO3_X1_TYPE0_NOERI)
      (sa, ia, W3c_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(w,a0) += (    1.00000000) D1(a1,a0) Fc1(w,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,a0,i,a) W4(w,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W4ca_sigma_ccov_coov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccov_coov_no4_x0_type0_noeri,G_IF_SIGMA_CCOV_COOV_NO4_X0_TYPE0_NOERI)
    (W4ca_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no4_x1_type0_noeri,G_IF_SIGMA_CCOV_COOV_NO4_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W4ca_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W5(x,a0) += (    1.00000000) D1(a1,a0) Fc1(x,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,a0,i,a) W5(x,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W5ca_sigma_ccov_coov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccov_coov_no5_x0_type0_noeri,G_IF_SIGMA_CCOV_COOV_NO5_X0_TYPE0_NOERI)
    (W5ca_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no5_x1_type0_noeri,G_IF_SIGMA_CCOV_COOV_NO5_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W5ca_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  } // End Femto
  //*-- FEMTO ends --//*

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
  // -- Title : sigma_ccov_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W7(c0,i,a2,a) += (    1.00000000) T2(c0,a1,a0,a) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) V2(c0,x,w,a2) W7(c0,i,a2,a) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W7aa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no0_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO0_X0_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, T2b.cptr(), W7aa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_coov_no0_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO0_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, V2_sym.cptr(), W7aa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W8(x,i,a1,a0) += (    1.00000000) V2(x,a3,a4,a2) D3(i,a0,a4,a2,a1,a3) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(w,a1,a0,a) W8(x,i,a1,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W8aaa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_coov_no1_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO1_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W8aaa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no1_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO1_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W8aaa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W9(w,i,a1,a0) += (    1.00000000) V2(w,a3,a4,a2) D3(i,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,a1,a0,a) W9(w,i,a1,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9aaa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_coov_no2_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO2_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W9aaa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no2_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO2_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W9aaa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W10(c0,a) += (    1.00000000) T2(c0,a1,a0,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) V2(c0,w,x,i) W10(c0,a) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    double W10_sigma_ccov_coov(0);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no3_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO3_X0_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, T2b.cptr(), &W10_sigma_ccov_coov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_coov_no3_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO3_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, V2_sym.cptr(), &W10_sigma_ccov_coov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W11(c0,a) += (    1.00000000) T2(c0,a1,a0,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) V2(c0,x,w,i) W11(c0,a) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    double W11_sigma_ccov_coov(0);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no4_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO4_X0_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, T2b.cptr(), &W11_sigma_ccov_coov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_coov_no4_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO4_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, V2_sym.cptr(), &W11_sigma_ccov_coov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W12(x,a1,a0,i) += (    1.00000000) V2(x,i,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(w,a1,a0,a) W12(x,a1,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W12aaa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_coov_no5_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO5_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W12aaa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no5_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO5_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W12aaa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W13(w,a1,a0,i) += (    1.00000000) V2(w,i,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,a1,a0,a) W13(w,a1,a0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W13aaa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_coov_no6_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO6_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W13aaa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no6_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO6_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W13aaa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W14(x,a1,a0,i) += (    1.00000000) V2(x,a3,i,a2) D2(a3,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(w,a0,a1,a) W14(x,a1,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W14aaa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_coov_no7_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO7_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W14aaa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no7_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO7_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W14aaa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W15(w,a1,a0,i) += (    1.00000000) V2(w,a3,i,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,a0,a1,a) W15(w,a1,a0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W15aaa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_coov_no8_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO8_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W15aaa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no8_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO8_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W15aaa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W16(x,w,c0,a0) += (    1.00000000) V2(x,c0,w,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) T2(c0,a0,i,a) W16(x,w,c0,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W16cca_sigma_ccov_coov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_coov_no9_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO9_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W16cca_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no9_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO9_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W16cca_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W17(x,w,c0,a0) += (    1.00000000) V2(x,a1,w,c0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(c0,a0,i,a) W17(x,w,c0,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W17cca_sigma_ccov_coov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_coov_no10_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO10_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W17cca_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no10_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO10_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W17cca_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W18(w,a0) += (    1.00000000) V2(w,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,a0,i,a) W18(w,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W18a_sigma_ccov_coov(orz::mr::sizeof_sympack_Xa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_coov_no11_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO11_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W18a_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no11_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO11_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W18a_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W19(x,a0) += (    1.00000000) V2(x,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,a0,i,a) W19(x,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W19a_sigma_ccov_coov(orz::mr::sizeof_sympack_Xa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_coov_no12_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO12_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W19a_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no12_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO12_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W19a_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W20(x,i,a0,a1) += (    1.00000000) V2(x,a2,a3,a1) D2(i,a3,a0,a2) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(w,a0,a1,a) W20(x,i,a0,a1) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W20aaa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_coov_no13_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO13_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W20aaa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no13_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO13_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W20aaa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W21(w,i,a0,a1) += (    1.00000000) V2(w,a2,a3,a1) D2(i,a2,a0,a3) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,a0,a1,a) W21(w,i,a0,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W21aaa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_coov_no14_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO14_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W21aaa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no14_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO14_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W21aaa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W22(x,a0,i,a1) += (    1.00000000) V2(x,i,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(w,a0,a1,a) W22(x,a0,i,a1) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W22aaa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_coov_no15_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO15_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W22aaa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no15_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO15_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W22aaa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W23(w,a0,i,a1) += (    1.00000000) V2(w,i,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,a0,a1,a) W23(w,a0,i,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W23aaa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_coov_no16_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO16_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W23aaa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no16_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO16_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W23aaa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W24(w,a0,i,a1) += (    1.00000000) V2(w,a2,i,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,a0,a1,a) W24(w,a0,i,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W24aaa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_coov_no17_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO17_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W24aaa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no17_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO17_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W24aaa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W25(x,a0,i,a1) += (    1.00000000) V2(x,a2,i,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,a0,a1,a) W25(x,a0,i,a1) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W25aaa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_coov_no18_x0_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO18_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W25aaa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no18_x1_type1_eri_c,G_IF_SIGMA_CCOV_COOV_NO18_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W25aaa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
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
  // -- Title : sigma_ccov_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W6(c0,i,a2,a) += (    1.00000000) T2(c0,a1,a0,a) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) V2(a2,x,w,c0) W6(c0,i,a2,a) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W6ca_sigma_ccov_coov(orz::mr::sizeof_sympack_Xca(symblockinfo, sa2^sa));
    FC_FUNC(g_if_sigma_ccov_coov_no0_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOV_NO0_X0_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, T2b.cptr(), W6ca_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no0_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOV_NO0_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W6ca_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
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
  // -- Title : sigma_ccov_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W26(w,i,a2,v0) += (    1.00000000) T2(w,a1,a0,v0) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) V2(v0,a2,x,a) W26(w,i,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W26caa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccov_coov_no0_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO0_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W26caa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no0_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO0_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W26caa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W27(x,i,a2,v0) += (    1.00000000) T2(x,a1,a0,v0) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) V2(v0,a2,w,a) W27(x,i,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W27caa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccov_coov_no1_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO1_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W27caa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no1_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO1_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W27caa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W28(w,i,a2,v0) += (    1.00000000) T2(w,a1,a0,v0) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) V2(v0,a,x,a2) W28(w,i,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W28caa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccov_coov_no2_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO2_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W28caa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no2_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO2_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W28caa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W29(x,i,a2,v0) += (    1.00000000) T2(x,a1,a0,v0) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) V2(v0,a,w,a2) W29(x,i,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W29caa_sigma_ccov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccov_coov_no3_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO3_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W29caa_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no3_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO3_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W29caa_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W30(w,v0) += (    1.00000000) T2(w,a1,a0,v0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) V2(v0,a,x,i) W30(w,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W30c_sigma_ccov_coov(orz::mr::sizeof_sympack_Xc(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccov_coov_no4_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO4_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W30c_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no4_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO4_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W30c_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W31(x,v0) += (    1.00000000) T2(x,a1,a0,v0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) V2(v0,a,w,i) W31(x,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W31c_sigma_ccov_coov(orz::mr::sizeof_sympack_Xc(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccov_coov_no5_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO5_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W31c_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no5_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO5_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W31c_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W32(x,a0,v0,a) += (    1.00000000) V2(a,x,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(w,a0,i,v0) W32(x,a0,v0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    orz::DTensor W32ca_sigma_ccov_coov(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_coov_no6_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO6_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W32ca_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_coov_no6_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO6_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W32ca_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W33(w,a0,v0,a) += (    1.00000000) V2(a,w,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,a0,i,v0) W33(w,a0,v0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    orz::DTensor W33ca_sigma_ccov_coov(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_coov_no7_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO7_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W33ca_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_coov_no7_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO7_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W33ca_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W34(w,a0,v0,a) += (    1.00000000) V2(v0,a,w,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,a0,i,v0) W34(w,a0,v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W34ca_sigma_ccov_coov(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_coov_no8_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO8_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W34ca_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_coov_no8_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO8_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W34ca_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W35(x,a0,v0,a) += (    1.00000000) V2(v0,a,x,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,a0,i,v0) W35(x,a0,v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W35ca_sigma_ccov_coov(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_coov_no9_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO9_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W35ca_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_coov_no9_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO9_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W35ca_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W36(w,v0) += (    1.00000000) T2(w,a1,a0,v0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) V2(v0,i,x,a) W36(w,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W36c_sigma_ccov_coov(orz::mr::sizeof_sympack_Xc(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccov_coov_no10_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO10_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W36c_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no10_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO10_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W36c_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W37(x,v0) += (    1.00000000) T2(x,a1,a0,v0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) V2(v0,i,w,a) W37(x,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W37c_sigma_ccov_coov(orz::mr::sizeof_sympack_Xc(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccov_coov_no11_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO11_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W37c_sigma_ccov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_coov_no11_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOV_NO11_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W37c_sigma_ccov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccov_coov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
