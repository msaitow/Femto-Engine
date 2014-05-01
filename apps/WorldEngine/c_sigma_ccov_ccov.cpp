                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ccov_ccov.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//      _/_/_/_/                            _/             
//     _/        _/_/    _/_/_/  _/_/    _/_/_/_/    _/_/  
//    _/_/_/  _/_/_/_/  _/    _/    _/    _/      _/    _/ 
//   _/      _/        _/    _/    _/    _/      _/    _/  
//  _/        _/_/_/  _/    _/    _/      _/_/    _/_/     

//                                   Generated date : Sun Apr 20 10:26:24 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccov_ccov(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| S2(w,x,i,a) += (    1.00000000) Fc0 T2(x,w,a0,a) D1(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no0_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO0_X0_TYPE0_NOERI)
      (sa, ia, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -2.00000000) Fc0 T2(w,x,a0,a) D1(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no1_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO1_X0_TYPE0_NOERI)
      (sa, ia, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    4.00000000) Fc0 T2(w,x,i,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no2_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO2_X0_TYPE0_NOERI)
      (sa, ia, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -2.00000000) Fc0 T2(x,w,i,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no3_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO3_X0_TYPE0_NOERI)
      (sa, ia, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W0(w,c0,i,a) += (    1.00000000) T2(c0,w,a0,a) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) Fc1(x,c0) W0(w,c0,i,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W0cca_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no4_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO4_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W0cca_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no4_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO4_X1_TYPE0_NOERI)
      (sa, ia, W0cca_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W1(x,c0,i,a) += (    1.00000000) T2(c0,x,a0,a) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) Fc1(w,c0) W1(x,c0,i,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W1cca_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no5_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO5_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W1cca_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no5_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO5_X1_TYPE0_NOERI)
      (sa, ia, W1cca_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W2(x,c0,i,a) += (    1.00000000) T2(x,c0,a0,a) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) Fc1(w,c0) W2(x,c0,i,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W2cca_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no6_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO6_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W2cca_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no6_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO6_X1_TYPE0_NOERI)
      (sa, ia, W2cca_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W3(w,c0,i,a) += (    1.00000000) T2(w,c0,a0,a) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) Fc1(x,c0) W3(w,c0,i,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W3cca_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no7_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO7_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W3cca_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no7_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO7_X1_TYPE0_NOERI)
      (sa, ia, W3cca_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W4(i,a0) += (    1.00000000) D2(i,a0,a2,a1) Fc1(a2,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,w,a0,a) W4(i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W4aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccov_ccov_no8_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO8_X0_TYPE0_NOERI)
    (W4aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no8_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO8_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W4aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W5(i,a0) += (    1.00000000) D2(i,a0,a2,a1) Fc1(a2,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,x,a0,a) W5(i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W5aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccov_ccov_no9_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO9_X0_TYPE0_NOERI)
    (W5aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no9_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO9_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W5aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W6(a0,i) += (    1.00000000) D1(a1,a0) Fc1(i,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,w,a0,a) W6(a0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W6aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccov_ccov_no10_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO10_X0_TYPE0_NOERI)
    (W6aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no10_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO10_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W6aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W7(a0,i) += (    1.00000000) D1(a1,a0) Fc1(i,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,x,a0,a) W7(a0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W7aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccov_ccov_no11_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO11_X0_TYPE0_NOERI)
    (W7aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no11_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO11_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W7aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) T2(w,c0,i,a) Fc1(x,c0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no12_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO12_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    2.00000000) T2(x,c0,i,a) Fc1(w,c0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no13_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO13_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    2.00000000) T2(c0,w,i,a) Fc1(x,c0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no14_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO14_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) T2(c0,x,i,a) Fc1(w,c0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no15_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO15_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W8() += (    1.00000000) D1(a1,a0) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(w,x,i,a) W8() 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  double W8_sigma_ccov_ccov(0);
  FC_FUNC(g_if_sigma_ccov_ccov_no16_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO16_X0_TYPE0_NOERI)
    (&W8_sigma_ccov_ccov, nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no16_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO16_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), &W8_sigma_ccov_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W9() += (    1.00000000) D1(a1,a0) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,w,i,a) W9() 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  double W9_sigma_ccov_ccov(0);
  FC_FUNC(g_if_sigma_ccov_ccov_no17_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO17_X0_TYPE0_NOERI)
    (&W9_sigma_ccov_ccov, nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no17_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO17_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), &W9_sigma_ccov_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W10(i,a0) += (    1.00000000) D1(i,a1) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,w,a0,a) W10(i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W10aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccov_ccov_no18_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO18_X0_TYPE0_NOERI)
    (W10aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no18_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO18_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W10aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W11(i,a0) += (    1.00000000) D1(i,a1) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,x,a0,a) W11(i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W11aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccov_ccov_no19_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO19_X0_TYPE0_NOERI)
    (W11aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no19_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO19_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W11aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    4.00000000) T2(w,x,a0,a) Fc1(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no20_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO20_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -2.00000000) T2(x,w,a0,a) Fc1(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no21_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO21_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W12(x,w,i,v0) += (    1.00000000) T2(x,w,a0,v0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) Fc1(v0,a) W12(x,w,i,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    orz::DTensor W12cca_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sv0));
    FC_FUNC(g_if_sigma_ccov_ccov_no22_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO22_X0_TYPE0_NOERI)
      (sv0, iv0, T2b.cptr(), W12cca_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_ccov_ccov_no22_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO22_X1_TYPE0_NOERI)
        (sa, ia, sv0, iv0, W12cca_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
      retval.acc_amp2(ia, S2b);
    } // End ia
    } // End sa
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W13(x,w,i,v0) += (    1.00000000) T2(x,w,v0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) Fc1(v0,a) W13(x,w,i,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W13ccav_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xccav(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_ccov_no23_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO23_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W13ccav_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no23_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO23_X1_TYPE0_NOERI)
      (sa, ia, W13ccav_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    4.00000000) T2(x,w,v0,i) Fc1(v0,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_ccov_ccov_no24_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO24_X0_TYPE0_NOERI)
        (sa, ia, si, ii, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
      retval.acc_amp2(ia, S2b);
    } // End ia
    } // End sa
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -2.00000000) T2(x,w,i,v0) Fc1(v0,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_ccov_ccov_no25_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOV_NO25_X0_TYPE0_NOERI)
        (sa, ia, sv0, iv0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
      retval.acc_amp2(ia, S2b);
    } // End ia
    } // End sa
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
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
  // -- Title : sigma_ccov_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W14(x,w,a0,a) += (    1.00000000) V2(x,c1,w,c0) T2(c1,c0,a0,a) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) D1(i,a0) W14(x,w,a0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W14ca_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sx^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no0_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO0_X0_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), V2_sym.cptr(), W14ca_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no0_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO0_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, W14ca_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W15(x,w,a0,a) += (    1.00000000) V2(x,c1,w,c0) T2(c0,c1,a0,a) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) D1(i,a0) W15(x,w,a0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W15ca_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sx^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no1_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO1_X0_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), V2_sym.cptr(), W15ca_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no1_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO1_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, W15ca_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W16(w,c0,i,a0) += (    1.00000000) V2(w,a2,c0,a1) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(c0,x,a0,a) W16(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W16caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no2_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO2_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W16caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no2_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO2_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W16caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W17(x,c0,i,a0) += (    1.00000000) V2(x,a2,c0,a1) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) T2(c0,w,a0,a) W17(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W17caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no3_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO3_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W17caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no3_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO3_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W17caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W18(w,c0,i,a0) += (    1.00000000) V2(w,a2,c0,a1) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) T2(x,c0,a0,a) W18(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W18caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no4_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO4_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W18caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no4_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO4_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W18caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W19(x,c0,i,a0) += (    1.00000000) V2(x,a2,c0,a1) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) T2(w,c0,a0,a) W19(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W19caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no5_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO5_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W19caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no5_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO5_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W19caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W20(w,c0,i,a0) += (    1.00000000) V2(w,c0,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(c0,x,a0,a) W20(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W20caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no6_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO6_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W20caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no6_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO6_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W20caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W21(x,c0,i,a0) += (    1.00000000) V2(x,c0,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) T2(c0,w,a0,a) W21(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W21caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no7_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO7_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W21caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no7_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO7_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W21caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W22(w,c0,i,a0) += (    1.00000000) V2(w,c0,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) T2(x,c0,a0,a) W22(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W22caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no8_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO8_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W22caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no8_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO8_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W22caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W23(x,c0,i,a0) += (    1.00000000) V2(x,c0,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(w,c0,a0,a) W23(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W23caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no9_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO9_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W23caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no9_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO9_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W23caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W26(w,c0,a0,i) += (    1.00000000) V2(w,i,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(c0,x,a0,a) W26(w,c0,a0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W26caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no10_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO10_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W26caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no10_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO10_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W26caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W27(x,c0,a0,i) += (    1.00000000) V2(x,i,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(c0,w,a0,a) W27(x,c0,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W27caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no11_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO11_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W27caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no11_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO11_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W27caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W28(w,c0,a0,i) += (    1.00000000) V2(w,i,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(x,c0,a0,a) W28(w,c0,a0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W28caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no12_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO12_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W28caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no12_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO12_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W28caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W29(x,c0,a0,i) += (    1.00000000) V2(x,i,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) T2(w,c0,a0,a) W29(x,c0,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W29caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no13_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO13_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W29caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no13_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO13_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W29caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W30(w,c0,a0,i) += (    1.00000000) V2(w,c0,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(c0,x,a0,a) W30(w,c0,a0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W30caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no14_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO14_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W30caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no14_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO14_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W30caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W31(x,c0,a0,i) += (    1.00000000) V2(x,c0,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) T2(c0,w,a0,a) W31(x,c0,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W31caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no15_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO15_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W31caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no15_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO15_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W31caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W32(w,c0,a0,i) += (    1.00000000) V2(w,c0,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) T2(x,c0,a0,a) W32(w,c0,a0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W32caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no16_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO16_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W32caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no16_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO16_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W32caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W33(x,c0,a0,i) += (    1.00000000) V2(x,c0,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(w,c0,a0,a) W33(x,c0,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W33caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no17_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO17_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W33caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no17_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO17_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W33caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    4.00000000) V2(x,c0,w,c1) T2(c1,c0,i,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no18_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO18_X0_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -2.00000000) V2(x,c1,w,c0) T2(c1,c0,i,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no19_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO19_X0_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W36(x,c0) += (    1.00000000) V2(x,a1,c0,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(w,c0,i,a) W36(x,c0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W36c_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no20_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO20_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W36c_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no20_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO20_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W36c_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W37(w,c0) += (    1.00000000) V2(w,a1,c0,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) T2(x,c0,i,a) W37(w,c0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W37c_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no21_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO21_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W37c_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no21_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO21_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W37c_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W38(x,c0) += (    1.00000000) V2(x,a1,c0,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) T2(c0,w,i,a) W38(x,c0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W38c_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no22_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO22_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W38c_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no22_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO22_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W38c_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W39(w,c0) += (    1.00000000) V2(w,a1,c0,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(c0,x,i,a) W39(w,c0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W39c_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no23_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO23_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W39c_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no23_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO23_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W39c_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W40(x,c0) += (    1.00000000) V2(x,c0,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(w,c0,i,a) W40(x,c0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W40c_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no24_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO24_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W40c_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no24_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO24_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W40c_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W41(w,c0) += (    1.00000000) V2(w,c0,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(x,c0,i,a) W41(w,c0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W41c_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no25_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO25_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W41c_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no25_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO25_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W41c_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W42(x,c0) += (    1.00000000) V2(x,c0,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(c0,w,i,a) W42(x,c0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W42c_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no26_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO26_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W42c_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no26_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO26_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W42c_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W43(w,c0) += (    1.00000000) V2(w,c0,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(c0,x,i,a) W43(w,c0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W43c_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no27_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO27_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W43c_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no27_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO27_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W43c_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W46(w,c0,i,a0) += (    1.00000000) V2(w,a1,c0,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(c0,x,a0,a) W46(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W46caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no28_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO28_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W46caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no28_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO28_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W46caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W47(x,c0,i,a0) += (    1.00000000) V2(x,a1,c0,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(c0,w,a0,a) W47(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W47caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no29_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO29_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W47caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no29_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO29_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W47caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W48(w,c0,i,a0) += (    1.00000000) V2(w,a1,c0,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(x,c0,a0,a) W48(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W48caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no30_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO30_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W48caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no30_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO30_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W48caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W49(x,c0,i,a0) += (    1.00000000) V2(x,a1,c0,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) T2(w,c0,a0,a) W49(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W49caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no31_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO31_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W49caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no31_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO31_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W49caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W50(w,c0,i,a0) += (    1.00000000) V2(w,c0,a1,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(c0,x,a0,a) W50(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W50caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no32_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO32_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W50caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no32_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO32_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W50caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W51(x,c0,i,a0) += (    1.00000000) V2(x,c0,a1,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) T2(c0,w,a0,a) W51(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W51caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no33_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO33_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W51caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no33_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO33_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W51caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W52(w,c0,i,a0) += (    1.00000000) V2(w,c0,a1,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) T2(x,c0,a0,a) W52(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W52caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccov_no34_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO34_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W52caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no34_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO34_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W52caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W53(x,c0,i,a0) += (    1.00000000) V2(x,c0,a1,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(w,c0,a0,a) W53(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W53caa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccov_no35_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO35_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W53caa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no35_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO35_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W53caa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    8.00000000) V2(w,i,c0,a0) T2(c0,x,a0,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no36_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO36_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) V2(x,i,c0,a0) T2(c0,w,a0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no37_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO37_X0_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) V2(w,i,c0,a0) T2(x,c0,a0,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no38_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO38_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    2.00000000) V2(x,i,c0,a0) T2(w,c0,a0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no39_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO39_X0_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) V2(x,c0,i,a0) T2(w,c0,a0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no40_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO40_X0_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    2.00000000) V2(w,c0,i,a0) T2(x,c0,a0,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no41_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO41_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    2.00000000) V2(x,c0,i,a0) T2(c0,w,a0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no42_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO42_X0_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) V2(w,c0,i,a0) T2(c0,x,a0,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no43_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCOV_NO43_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@type(2).declaration(begin)
  // --  Title : sigma_ccov_ccov
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W24aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  orz::DTensor W25aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  double W44_sigma_ccov_ccov(0);
  double W45_sigma_ccov_ccov(0);
//-@type(2).declaration(end)

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
  // -- Title : sigma_ccov_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W24(i,a0) += (    1.00000000) V2(a4,a2,a3,a1) D3(i,a0,a4,a2,a3,a1) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_ccov_ccov_no0_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO0_X0_TYPE1_ERI_O)
    (sa4, ia4, V2_sym.cptr(), W24aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W25(i,a0) += (    1.00000000) V2(a4,a2,a3,a1) D3(i,a0,a4,a2,a3,a1) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_ccov_ccov_no1_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO1_X0_TYPE1_ERI_O)
    (sa4, ia4, V2_sym.cptr(), W25aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W34(a0,i) += (    1.00000000) V2(i,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,x,a0,a) W34(a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W34a_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, si));
  FC_FUNC(g_if_sigma_ccov_ccov_no2_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO2_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W34a_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no2_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO2_X1_TYPE1_ERI_O)
      (sa, ia, si, ii, T2b.cptr(), W34a_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W35(a0,i) += (    1.00000000) V2(i,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,w,a0,a) W35(a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W35a_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, si));
  FC_FUNC(g_if_sigma_ccov_ccov_no3_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO3_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W35a_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no3_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO3_X1_TYPE1_ERI_O)
      (sa, ia, si, ii, T2b.cptr(), W35a_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W44() += (    1.00000000) V2(a3,a1,a2,a0) D2(a3,a1,a2,a0) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_ccov_ccov_no4_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO4_X0_TYPE1_ERI_O)
    (sa3, ia3, V2_sym.cptr(), &W44_sigma_ccov_ccov, nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W45() += (    1.00000000) V2(a3,a1,a2,a0) D2(a3,a1,a2,a0) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_ccov_ccov_no5_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO5_X0_TYPE1_ERI_O)
    (sa3, ia3, V2_sym.cptr(), &W45_sigma_ccov_ccov, nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W54(i,a0) += (    1.00000000) V2(a0,a2,a3,a1) D2(i,a2,a3,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,w,a,a0) W54(i,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia0);
  orz::DTensor W54a_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sa0));
  FC_FUNC(g_if_sigma_ccov_ccov_no6_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO6_X0_TYPE1_ERI_O)
    (sa0, ia0, V2_sym.cptr(), W54a_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no6_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO6_X1_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, T2b.cptr(), W54a_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W55(i,a0) += (    1.00000000) V2(a0,a2,a3,a1) D2(i,a2,a3,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,w,a0,a) W55(i,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W55a_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sa0));
  FC_FUNC(g_if_sigma_ccov_ccov_no7_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO7_X0_TYPE1_ERI_O)
    (sa0, ia0, V2_sym.cptr(), W55a_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no7_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO7_X1_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, T2b.cptr(), W55a_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W56(i,a0) += (    1.00000000) V2(i,a1,a2,a0) D1(a2,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,x,a0,a) W56(i,a0) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W56a_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, si));
  FC_FUNC(g_if_sigma_ccov_ccov_no8_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO8_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W56a_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no8_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO8_X1_TYPE1_ERI_O)
      (sa, ia, si, ii, T2b.cptr(), W56a_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W57(i,a0) += (    1.00000000) V2(i,a1,a2,a0) D1(a2,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,w,a0,a) W57(i,a0) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W57a_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, si));
  FC_FUNC(g_if_sigma_ccov_ccov_no9_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO9_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W57a_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no9_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO9_X1_TYPE1_ERI_O)
      (sa, ia, si, ii, T2b.cptr(), W57a_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W58(i,a1) += (    1.00000000) V2(i,a1,a2,a0) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(w,x,a1,a) W58(i,a1) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W58a_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, si));
  FC_FUNC(g_if_sigma_ccov_ccov_no10_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO10_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W58a_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no10_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO10_X1_TYPE1_ERI_O)
      (sa, ia, si, ii, T2b.cptr(), W58a_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W59(i,a1) += (    1.00000000) V2(i,a1,a2,a0) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,w,a1,a) W59(i,a1) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W59a_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, si));
  FC_FUNC(g_if_sigma_ccov_ccov_no11_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO11_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W59a_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no11_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOV_NO11_X1_TYPE1_ERI_O)
      (sa, ia, si, ii, T2b.cptr(), W59a_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@type(2).contraction(begin)
  // -- Title : sigma_ccov_ccov
  //*-- Entering to take the type 2 contractions --*//

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -1.00000000) T2(w,x,a0,a) W24(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no0_x0_type2_eri_o,G_IF_SIGMA_CCOV_CCOV_NO0_X0_TYPE2_ERI_O)
      (sa, ia, T2b.cptr(), W24aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    0.50000000) T2(x,w,a0,a) W25(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no1_x0_type2_eri_o,G_IF_SIGMA_CCOV_CCOV_NO1_X0_TYPE2_ERI_O)
      (sa, ia, T2b.cptr(), W25aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    2.00000000) T2(w,x,i,a) W44() 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no2_x0_type2_eri_o,G_IF_SIGMA_CCOV_CCOV_NO2_X0_TYPE2_ERI_O)
      (sa, ia, T2b.cptr(), &W44_sigma_ccov_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -1.00000000) T2(x,w,i,a) W45() 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no3_x0_type2_eri_o,G_IF_SIGMA_CCOV_CCOV_NO3_X0_TYPE2_ERI_O)
      (sa, ia, T2b.cptr(), &W45_sigma_ccov_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
//-@type(2).contraction(end)
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
  // -- Title : sigma_ccov_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W60(w,x,a0,a) += (    1.00000000) V2(a,x,v0,c0) T2(w,c0,v0,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) D1(i,a0) W60(w,x,a0,a) 
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
    orz::DTensor W60cc_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no0_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO0_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W60cc_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no0_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO0_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W60cc_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W61(x,w,a0,a) += (    1.00000000) V2(a,w,v0,c0) T2(x,c0,v0,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) D1(i,a0) W61(x,w,a0,a) 
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
    orz::DTensor W61cc_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no1_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO1_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W61cc_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no1_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO1_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W61cc_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W62(x,w,a0,a) += (    1.00000000) V2(a,w,v0,c0) T2(x,c0,a0,v0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) D1(i,a0) W62(x,w,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W62cca_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccov_ccov_no2_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO2_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W62cca_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ccov_ccov_no2_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO2_X1_TYPE1_ERI_V)
    (sa, ia, W62cca_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W63(w,x,a0,a) += (    1.00000000) V2(a,x,v0,c0) T2(w,c0,a0,v0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) D1(i,a0) W63(w,x,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W63cca_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccov_ccov_no3_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO3_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W63cca_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ccov_ccov_no3_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO3_X1_TYPE1_ERI_V)
    (sa, ia, W63cca_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W64(i,a0,v0,a) += (    1.00000000) V2(v0,a1,a2,a) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(w,x,a0,v0) W64(i,a0,v0,a) 
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
    orz::DTensor W64aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no4_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO4_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W64aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no4_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO4_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W64aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W65(i,a0,v0,a) += (    1.00000000) V2(v0,a1,a2,a) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,w,a0,v0) W65(i,a0,v0,a) 
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
    orz::DTensor W65aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no5_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO5_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W65aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no5_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO5_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W65aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W66(x,w,a0,a) += (    1.00000000) V2(a,v0,w,c0) T2(x,c0,v0,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) D1(i,a0) W66(x,w,a0,a) 
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
    orz::DTensor W66cc_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no6_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO6_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W66cc_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no6_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO6_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W66cc_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W67(w,x,a0,a) += (    1.00000000) V2(a,v0,x,c0) T2(w,c0,v0,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) D1(i,a0) W67(w,x,a0,a) 
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
    orz::DTensor W67cc_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no7_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO7_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W67cc_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no7_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO7_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W67cc_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W68(x,w,a0,a) += (    1.00000000) V2(a,v0,w,c0) T2(x,c0,a0,v0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -1.00000000) D1(i,a0) W68(x,w,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W68cca_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccov_ccov_no8_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO8_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W68cca_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ccov_ccov_no8_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO8_X1_TYPE1_ERI_V)
    (sa, ia, W68cca_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W69(w,x,a0,a) += (    1.00000000) V2(a,v0,x,c0) T2(w,c0,a0,v0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) D1(i,a0) W69(w,x,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W69cca_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccov_ccov_no9_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO9_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W69cca_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ccov_ccov_no9_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO9_X1_TYPE1_ERI_V)
    (sa, ia, W69cca_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W70(i,a0,v0,a) += (    1.00000000) V2(v0,a,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,x,a0,v0) W70(i,a0,v0,a) 
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
    orz::DTensor W70aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no10_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO10_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W70aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no10_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO10_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W70aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W71(i,a0,v0,a) += (    1.00000000) V2(v0,a,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,w,a0,v0) W71(i,a0,v0,a) 
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
    orz::DTensor W71aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no11_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO11_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W71aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no11_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO11_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W71aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W72(a0,i,v0,a) += (    1.00000000) V2(v0,a,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,x,a0,v0) W72(a0,i,v0,a) 
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
    orz::DTensor W72aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no12_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO12_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W72aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no12_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO12_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W72aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W73(a0,i,v0,a) += (    1.00000000) V2(v0,a,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,w,a0,v0) W73(a0,i,v0,a) 
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
    orz::DTensor W73aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no13_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO13_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W73aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no13_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO13_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W73aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    8.00000000) V2(a,x,v0,c0) T2(w,c0,i,v0) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccov_ccov_no14_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO14_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) V2(a,w,v0,c0) T2(x,c0,i,v0) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccov_ccov_no15_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO15_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    2.00000000) V2(a,w,v0,c0) T2(x,c0,v0,i) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_ccov_no16_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO16_X0_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) V2(a,x,v0,c0) T2(w,c0,v0,i) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_ccov_no17_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO17_X0_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W74(v0,a) += (    1.00000000) V2(v0,a0,a1,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,w,i,v0) W74(v0,a) 
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
    double W74_sigma_ccov_ccov(0);
    FC_FUNC(g_if_sigma_ccov_ccov_no18_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO18_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), &W74_sigma_ccov_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no18_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO18_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), &W74_sigma_ccov_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W75(v0,a) += (    1.00000000) V2(v0,a0,a1,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,x,i,v0) W75(v0,a) 
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
    double W75_sigma_ccov_ccov(0);
    FC_FUNC(g_if_sigma_ccov_ccov_no19_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO19_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), &W75_sigma_ccov_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no19_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO19_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), &W75_sigma_ccov_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) V2(v0,a,x,c0) T2(w,c0,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no20_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO20_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    2.00000000) V2(v0,a,w,c0) T2(x,c0,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no21_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO21_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    2.00000000) V2(v0,a,x,c0) T2(c0,w,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no22_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO22_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) V2(v0,a,w,c0) T2(c0,x,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no23_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO23_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W76(v0,a) += (    1.00000000) V2(v0,a,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(w,x,i,v0) W76(v0,a) 
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
    double W76_sigma_ccov_ccov(0);
    FC_FUNC(g_if_sigma_ccov_ccov_no24_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO24_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), &W76_sigma_ccov_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no24_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO24_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), &W76_sigma_ccov_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W77(v0,a) += (    1.00000000) V2(v0,a,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,w,i,v0) W77(v0,a) 
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
    double W77_sigma_ccov_ccov(0);
    FC_FUNC(g_if_sigma_ccov_ccov_no25_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO25_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), &W77_sigma_ccov_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no25_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO25_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), &W77_sigma_ccov_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W78(a0,i,v0,a) += (    1.00000000) V2(v0,i,a1,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(w,x,a0,v0) W78(a0,i,v0,a) 
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
    orz::DTensor W78aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no26_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO26_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W78aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no26_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO26_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W78aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W79(a0,i,v0,a) += (    1.00000000) V2(v0,i,a1,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,w,a0,v0) W79(a0,i,v0,a) 
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
    orz::DTensor W79aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no27_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO27_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W79aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no27_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO27_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W79aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W80(i,a1,v0,a) += (    1.00000000) V2(v0,a0,a1,a) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,w,a1,v0) W80(i,a1,v0,a) 
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
    orz::DTensor W80aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no28_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO28_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W80aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no28_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO28_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W80aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W81(i,a1,v0,a) += (    1.00000000) V2(v0,a0,a1,a) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(w,x,a1,v0) W81(i,a1,v0,a) 
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
    orz::DTensor W81aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no29_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO29_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W81aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no29_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO29_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W81aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W82(i,a0,v0,a) += (    1.00000000) V2(v0,a,a1,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,x,a0,v0) W82(i,a0,v0,a) 
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
    orz::DTensor W82aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no30_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO30_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W82aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no30_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO30_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W82aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W83(i,a0,v0,a) += (    1.00000000) V2(v0,a,a1,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,w,a0,v0) W83(i,a0,v0,a) 
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
    orz::DTensor W83aa_sigma_ccov_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ccov_ccov_no31_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO31_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W83aa_sigma_ccov_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccov_no31_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO31_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W83aa_sigma_ccov_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    4.00000000) V2(v0,a,i,a0) T2(w,x,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no32_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO32_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -2.00000000) V2(v0,a,i,a0) T2(x,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no33_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO33_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    4.00000000) V2(v0,i,a0,a) T2(x,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no34_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO34_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -2.00000000) V2(v0,i,a0,a) T2(w,x,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccov_no35_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOV_NO35_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccov_ccov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
