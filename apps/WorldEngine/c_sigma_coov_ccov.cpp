                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_coov_ccov.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
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

//                                   Generated date : Sun Apr 20 10:26:19 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_coov_ccov(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(w,a1,a0,a) += (    1.00000000) T2(c0,w,a0,a) Fc1(c0,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D2(j,i,a1,a0) W0(w,a1,a0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W0caa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no0_x0_type0_noeri,G_IF_SIGMA_COOV_CCOV_NO0_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W0caa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_ccov_no0_x1_type0_noeri,G_IF_SIGMA_COOV_CCOV_NO0_X1_TYPE0_NOERI)
      (sa, ia, W0caa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(w,a1,a0,a) += (    1.00000000) T2(w,c0,a0,a) Fc1(c0,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D2(j,a0,a1,i) W1(w,a1,a0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W1caa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no1_x0_type0_noeri,G_IF_SIGMA_COOV_CCOV_NO1_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W1caa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_ccov_no1_x1_type0_noeri,G_IF_SIGMA_COOV_CCOV_NO1_X1_TYPE0_NOERI)
      (sa, ia, W1caa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(c0,i) += (    1.00000000) D1(i,a0) Fc1(c0,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) T2(c0,w,j,a) W2(c0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W2ca_sigma_coov_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_coov_ccov_no2_x0_type0_noeri,G_IF_SIGMA_COOV_CCOV_NO2_X0_TYPE0_NOERI)
    (W2ca_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no2_x1_type0_noeri,G_IF_SIGMA_COOV_CCOV_NO2_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W2ca_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(c0,i) += (    1.00000000) D1(i,a0) Fc1(c0,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) T2(w,c0,j,a) W3(c0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W3ca_sigma_coov_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_coov_ccov_no3_x0_type0_noeri,G_IF_SIGMA_COOV_CCOV_NO3_X0_TYPE0_NOERI)
    (W3ca_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no3_x1_type0_noeri,G_IF_SIGMA_COOV_CCOV_NO3_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W3ca_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(w,a) += (    1.00000000) T2(c0,w,a0,a) Fc1(c0,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) D1(j,i) W4(w,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W4c_sigma_coov_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no4_x0_type0_noeri,G_IF_SIGMA_COOV_CCOV_NO4_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W4c_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_ccov_no4_x1_type0_noeri,G_IF_SIGMA_COOV_CCOV_NO4_X1_TYPE0_NOERI)
      (sa, ia, W4c_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W5(w,a) += (    1.00000000) T2(w,c0,a0,a) Fc1(c0,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D1(j,i) W5(w,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W5c_sigma_coov_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no5_x0_type0_noeri,G_IF_SIGMA_COOV_CCOV_NO5_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W5c_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_ccov_no5_x1_type0_noeri,G_IF_SIGMA_COOV_CCOV_NO5_X1_TYPE0_NOERI)
      (sa, ia, W5c_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_coov_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W6(w,a0,a1,a) += (    1.00000000) V2(w,c0,c1,a1) T2(c1,c0,a0,a) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) D2(j,i,a1,a0) W6(w,a0,a1,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W6aa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no0_x0_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO0_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), W6aa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_ccov_no0_x1_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO0_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, W6aa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W7(w,a0,a1,a) += (    1.00000000) V2(w,c0,c1,a1) T2(c0,c1,a0,a) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) D2(j,a0,a1,i) W7(w,a0,a1,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W7aa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no1_x0_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO1_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), W7aa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_ccov_no1_x1_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO1_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, W7aa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W8(c0,j,i,a0) += (    1.00000000) V2(c0,a2,a3,a1) D3(j,i,a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) T2(c0,w,a0,a) W8(c0,j,i,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W8aaa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_coov_ccov_no2_x0_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO2_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W8aaa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no2_x1_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO2_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, T2b.cptr(), W8aaa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W9(c0,j,a0,i) += (    1.00000000) V2(c0,a2,a3,a1) D3(j,a0,a3,a1,a2,i) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) T2(w,c0,a0,a) W9(c0,j,a0,i) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9aaa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_coov_ccov_no3_x0_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO3_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W9aaa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no3_x1_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO3_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, T2b.cptr(), W9aaa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W12(w,c1,c0,i) += (    1.00000000) V2(w,c0,c1,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) T2(c1,c0,j,a) W12(w,c1,c0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W12cca_sigma_coov_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sw));
  FC_FUNC(g_if_sigma_coov_ccov_no4_x0_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO4_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W12cca_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no4_x1_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO4_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W12cca_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W13(w,c1,c0,i) += (    1.00000000) V2(w,c0,c1,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    2.00000000) T2(c0,c1,j,a) W13(w,c1,c0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W13cca_sigma_coov_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sw));
  FC_FUNC(g_if_sigma_coov_ccov_no5_x0_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO5_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W13cca_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no5_x1_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO5_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W13cca_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W14(c0,i) += (    1.00000000) V2(c0,a1,a2,a0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) T2(w,c0,j,a) W14(c0,i) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W14a_sigma_coov_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_coov_ccov_no6_x0_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO6_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W14a_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no6_x1_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO6_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, T2b.cptr(), W14a_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W15(c0,i) += (    1.00000000) V2(c0,a1,a2,a0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) T2(c0,w,j,a) W15(c0,i) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W15a_sigma_coov_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_coov_ccov_no7_x0_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO7_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W15a_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no7_x1_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO7_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, T2b.cptr(), W15a_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W16(w,a) += (    1.00000000) V2(w,c0,c1,a0) T2(c1,c0,a0,a) 
  // |-- [    1] --| S2(w,i,j,a) += (    2.00000000) D1(j,i) W16(w,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    double W16_sigma_coov_ccov(0);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no8_x0_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO8_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), &W16_sigma_coov_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_ccov_no8_x1_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO8_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, &W16_sigma_coov_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W17(w,a) += (    1.00000000) V2(w,c0,c1,a0) T2(c0,c1,a0,a) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) D1(j,i) W17(w,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    double W17_sigma_coov_ccov(0);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no9_x0_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO9_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), &W17_sigma_coov_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_ccov_no9_x1_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO9_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, &W17_sigma_coov_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W18(c0,j,i,a0) += (    1.00000000) V2(c0,a1,a2,a0) D2(j,i,a1,a2) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) T2(c0,w,a0,a) W18(c0,j,i,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W18aaa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_coov_ccov_no10_x0_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO10_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W18aaa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no10_x1_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO10_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, T2b.cptr(), W18aaa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W19(c0,j,i,a0) += (    1.00000000) V2(c0,a1,a2,a0) D2(j,a2,a1,i) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) T2(w,c0,a0,a) W19(c0,j,i,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W19aaa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_coov_ccov_no11_x0_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO11_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W19aaa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no11_x1_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO11_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, T2b.cptr(), W19aaa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W20(c0,j,i,a1) += (    1.00000000) V2(c0,a1,a2,a0) D2(j,i,a2,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) T2(c0,w,a1,a) W20(c0,j,i,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W20aaa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_coov_ccov_no12_x0_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO12_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W20aaa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no12_x1_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO12_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, T2b.cptr(), W20aaa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W21(c0,j,i,a1) += (    1.00000000) V2(c0,a1,a2,a0) D2(j,i,a2,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) T2(w,c0,a1,a) W21(c0,j,i,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W21aaa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_coov_ccov_no13_x0_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO13_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W21aaa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no13_x1_type1_eri_c,G_IF_SIGMA_COOV_CCOV_NO13_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, T2b.cptr(), W21aaa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_coov_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W10(c0,i,a0,j) += (    1.00000000) V2(j,a2,c0,a1) D2(i,a2,a0,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) T2(c0,w,a0,a) W10(c0,i,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W10caa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_coov_ccov_no0_x0_type1_eri_o,G_IF_SIGMA_COOV_CCOV_NO0_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W10caa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no0_x1_type1_eri_o,G_IF_SIGMA_COOV_CCOV_NO0_X1_TYPE1_ERI_O)
      (sa, ia, sj, ij, T2b.cptr(), W10caa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W11(c0,i,a0,j) += (    1.00000000) V2(j,a2,c0,a1) D2(i,a1,a0,a2) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) T2(w,c0,a0,a) W11(c0,i,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W11caa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_coov_ccov_no1_x0_type1_eri_o,G_IF_SIGMA_COOV_CCOV_NO1_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W11caa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no1_x1_type1_eri_o,G_IF_SIGMA_COOV_CCOV_NO1_X1_TYPE1_ERI_O)
      (sa, ia, sj, ij, T2b.cptr(), W11caa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W22(c0,i,j,a0) += (    1.00000000) V2(j,a1,c0,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) T2(c0,w,a0,a) W22(c0,i,j,a0) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W22caa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_coov_ccov_no2_x0_type1_eri_o,G_IF_SIGMA_COOV_CCOV_NO2_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W22caa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no2_x1_type1_eri_o,G_IF_SIGMA_COOV_CCOV_NO2_X1_TYPE1_ERI_O)
      (sa, ia, sj, ij, T2b.cptr(), W22caa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W23(c0,i,j,a0) += (    1.00000000) V2(j,a1,c0,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) T2(w,c0,a0,a) W23(c0,i,j,a0) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W23caa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_coov_ccov_no3_x0_type1_eri_o,G_IF_SIGMA_COOV_CCOV_NO3_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W23caa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no3_x1_type1_eri_o,G_IF_SIGMA_COOV_CCOV_NO3_X1_TYPE1_ERI_O)
      (sa, ia, sj, ij, T2b.cptr(), W23caa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W24(c0,i,j,a1) += (    1.00000000) V2(j,a1,c0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) T2(w,c0,a1,a) W24(c0,i,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W24caa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_coov_ccov_no4_x0_type1_eri_o,G_IF_SIGMA_COOV_CCOV_NO4_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W24caa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no4_x1_type1_eri_o,G_IF_SIGMA_COOV_CCOV_NO4_X1_TYPE1_ERI_O)
      (sa, ia, sj, ij, T2b.cptr(), W24caa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W25(c0,i,j,a1) += (    1.00000000) V2(j,a1,c0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) T2(c0,w,a1,a) W25(c0,i,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W25caa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_coov_ccov_no5_x0_type1_eri_o,G_IF_SIGMA_COOV_CCOV_NO5_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W25caa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_ccov_no5_x1_type1_eri_o,G_IF_SIGMA_COOV_CCOV_NO5_X1_TYPE1_ERI_O)
      (sa, ia, sj, ij, T2b.cptr(), W25caa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_coov_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W26(w,a0,a1,a) += (    1.00000000) V2(a,a1,v0,c0) T2(w,c0,v0,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D2(j,a0,a1,i) W26(w,a0,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W26ca_sigma_coov_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_coov_ccov_no0_x0_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO0_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W26ca_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_ccov_no0_x1_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO0_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W26ca_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W27(w,a0,a1,a) += (    1.00000000) V2(a,a1,v0,c0) T2(w,c0,a0,v0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) D2(j,a0,a1,i) W27(w,a0,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W27caa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_coov_ccov_no1_x0_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO1_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W27caa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_coov_ccov_no1_x1_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO1_X1_TYPE1_ERI_V)
    (sa, ia, W27caa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W28(w,a0,a1,a) += (    1.00000000) V2(a,v0,c0,a1) T2(w,c0,v0,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D2(j,i,a1,a0) W28(w,a0,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W28ca_sigma_coov_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_coov_ccov_no2_x0_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO2_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W28ca_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_ccov_no2_x1_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO2_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W28ca_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W29(w,a0,a1,a) += (    1.00000000) V2(a,v0,c0,a1) T2(w,c0,a0,v0) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D2(j,a0,a1,i) W29(w,a0,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W29caa_sigma_coov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_coov_ccov_no3_x0_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO3_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W29caa_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_coov_ccov_no3_x1_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO3_X1_TYPE1_ERI_V)
    (sa, ia, W29caa_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W30(c0,i,v0,a) += (    1.00000000) V2(v0,c0,a0,a) D1(i,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    4.00000000) T2(w,c0,j,v0) W30(c0,i,v0,a) 
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
    orz::DTensor W30ca_sigma_coov_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_coov_ccov_no4_x0_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO4_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W30ca_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_ccov_no4_x1_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO4_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W30ca_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W31(c0,i,v0,a) += (    1.00000000) V2(v0,c0,a0,a) D1(i,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) T2(c0,w,j,v0) W31(c0,i,v0,a) 
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
    orz::DTensor W31ca_sigma_coov_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_coov_ccov_no5_x0_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO5_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W31ca_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_ccov_no5_x1_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO5_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W31ca_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W32(c0,i,v0,a) += (    1.00000000) V2(v0,a,c0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) T2(w,c0,j,v0) W32(c0,i,v0,a) 
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
    orz::DTensor W32ca_sigma_coov_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_coov_ccov_no6_x0_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO6_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W32ca_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_ccov_no6_x1_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO6_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W32ca_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W33(c0,i,v0,a) += (    1.00000000) V2(v0,a,c0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) T2(c0,w,j,v0) W33(c0,i,v0,a) 
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
    orz::DTensor W33ca_sigma_coov_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_coov_ccov_no7_x0_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO7_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W33ca_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_ccov_no7_x1_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO7_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W33ca_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W34(w,a) += (    1.00000000) V2(a,a0,v0,c0) T2(w,c0,a0,v0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) D1(j,i) W34(w,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W34c_sigma_coov_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_coov_ccov_no8_x0_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO8_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W34c_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_coov_ccov_no8_x1_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO8_X1_TYPE1_ERI_V)
    (sa, ia, W34c_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W35(w,a) += (    1.00000000) V2(a,a0,v0,c0) T2(w,c0,v0,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D1(j,i) W35(w,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W35c_sigma_coov_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sa));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_coov_ccov_no9_x0_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO9_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W35c_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  FC_FUNC(g_if_sigma_coov_ccov_no9_x1_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO9_X1_TYPE1_ERI_V)
    (sa, ia, W35c_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W36(w,a) += (    1.00000000) V2(a,v0,c0,a0) T2(w,c0,v0,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) D1(j,i) W36(w,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W36c_sigma_coov_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sa));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_coov_ccov_no10_x0_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO10_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W36c_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  FC_FUNC(g_if_sigma_coov_ccov_no10_x1_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO10_X1_TYPE1_ERI_V)
    (sa, ia, W36c_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W37(w,a) += (    1.00000000) V2(a,v0,c0,a0) T2(w,c0,a0,v0) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D1(j,i) W37(w,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W37c_sigma_coov_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_coov_ccov_no11_x0_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO11_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W37c_sigma_coov_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_coov_ccov_no11_x1_type1_eri_v,G_IF_SIGMA_COOV_CCOV_NO11_X1_TYPE1_ERI_V)
    (sa, ia, W37c_sigma_coov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_coov_ccov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
