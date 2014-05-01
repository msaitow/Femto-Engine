                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ocov_coov.h>                                  
                                                                                
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

//                                   Generated date : Sun Apr 20 10:26:11 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ocov_coov(const orz::mr::Input &ctinp,                                    
                                  const orz::mr::SymBlockInfo &symblockinfo,                                 
                                  const orz::mr::HintMO &hintmo,                                             
                                  const int alloc_type,                                                      
                                  const orz::mr::RdmPack &rdmPack,                                           
                                  const orz::DTensor &rdm4,                                                  
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
  orz::DTensor rdm4_sym;                                                                                         
  orz::DTensor rdm4_ij_sliced(ctinp.use_d4cum_of() ? nocc*nocc*nocc*nocc*nocc*nocc : 0);                         
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
  // |-- [    0] --| S2(i,w,j,a) += (   -1.00000000) Fc0 T2(w,a1,a0,a) D2(j,i,a1,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no0_x0_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO0_X0_TYPE0_NOERI)
      (sa, ia, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(i,w,j,a) += (   -1.00000000) Fc0 T2(w,a0,j,a) D1(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no1_x0_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO1_X0_TYPE0_NOERI)
      (sa, ia, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W0(c0,j,i,a) += (    1.00000000) T2(c0,a1,a0,a) D2(j,i,a1,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (    1.00000000) Fc1(w,c0) W0(c0,j,i,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W0caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no2_x0_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO2_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W0caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_coov_no2_x1_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO2_X1_TYPE0_NOERI)
      (sa, ia, W0caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W1(j,a1,i,a0) += (    1.00000000) D3(j,i,a3,a2,a1,a0) Fc1(a3,a2) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) T2(w,a1,a0,a) W1(j,a1,i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1aaaa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ocov_coov_no3_x0_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO3_X0_TYPE0_NOERI)
    (W1aaaa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no3_x1_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO3_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W1aaaa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W2(i,a1,a0,j) += (    1.00000000) D2(i,a2,a1,a0) Fc1(j,a2) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) T2(w,a0,a1,a) W2(i,a1,a0,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W2aaaa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ocov_coov_no4_x0_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO4_X0_TYPE0_NOERI)
    (W2aaaa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no4_x1_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO4_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W2aaaa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W3(c0,i,j,a) += (    1.00000000) T2(c0,a0,j,a) D1(i,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (    1.00000000) Fc1(w,c0) W3(c0,i,j,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W3caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no5_x0_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO5_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W3caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_coov_no5_x1_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO5_X1_TYPE0_NOERI)
      (sa, ia, W3caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W4(i,a0) += (    1.00000000) D2(i,a0,a2,a1) Fc1(a2,a1) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) T2(w,a0,j,a) W4(i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W4aa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ocov_coov_no6_x0_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO6_X0_TYPE0_NOERI)
    (W4aa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no6_x1_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO6_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W4aa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W5(j,a0,i,a1) += (    1.00000000) D2(j,i,a0,a2) Fc1(a2,a1) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) T2(w,a0,a1,a) W5(j,a0,i,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W5aaaa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ocov_coov_no7_x0_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO7_X0_TYPE0_NOERI)
    (W5aaaa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no7_x1_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO7_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W5aaaa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W6(w,i,a1,a) += (    1.00000000) T2(w,a0,a1,a) D1(i,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) Fc1(j,a1) W6(w,i,a1,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W6caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no8_x0_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO8_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W6caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_coov_no8_x1_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO8_X1_TYPE0_NOERI)
      (sa, ia, W6caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W7(w,j,i,v0) += (    1.00000000) T2(w,a1,a0,v0) D2(j,i,a1,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) Fc1(v0,a) W7(w,j,i,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    orz::DTensor W7caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
    FC_FUNC(g_if_sigma_ocov_coov_no9_x0_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO9_X0_TYPE0_NOERI)
      (sv0, iv0, T2b.cptr(), W7caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_ocov_coov_no9_x1_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO9_X1_TYPE0_NOERI)
        (sa, ia, sv0, iv0, W7caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   10] -- 
  // |-- [    0] --| W8(w,i,j,v0) += (    1.00000000) T2(w,a0,j,v0) D1(i,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) Fc1(v0,a) W8(w,i,j,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    orz::DTensor W8caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
    FC_FUNC(g_if_sigma_ocov_coov_no10_x0_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO10_X0_TYPE0_NOERI)
      (sv0, iv0, T2b.cptr(), W8caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_ocov_coov_no10_x1_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO10_X1_TYPE0_NOERI)
        (sa, ia, sv0, iv0, W8caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   11] -- 
  // |-- [    0] --| S2(i,w,j,a) += (   -0.50000000) T2(w,a1,a0,a) C4(a0,a1,i,j) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no11_x0_type0_noeri,G_IF_SIGMA_OCOV_COOV_NO11_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ocov_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W9(w,c0,j,a1,i,a0) += (    1.00000000) V2(w,a3,c0,a2) D3(j,i,a2,a0,a1,a3) 
  // |-- [    1] --| S2(i,w,j,a) += (    1.00000000) T2(c0,a1,a0,a) W9(w,c0,j,a1,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9caaaa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ocov_coov_no0_x0_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO0_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W9caaaa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no0_x1_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO0_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W9caaaa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W10(w,c0,j,a1,i,a0) += (    1.00000000) V2(w,c0,a3,a2) D3(j,i,a3,a2,a1,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (    1.00000000) T2(c0,a1,a0,a) W10(w,c0,j,a1,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W10caaaa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ocov_coov_no1_x0_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO1_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W10caaaa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no1_x1_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO1_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W10caaaa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W11(c0,i,a2,a) += (    1.00000000) T2(c0,a0,a1,a) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(i,w,j,a) += (    1.00000000) V2(c0,a2,w,j) W11(c0,i,a2,a) 
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
    orz::DTensor W11aa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no2_x0_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO2_X0_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, T2b.cptr(), W11aa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_coov_no2_x1_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO2_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, V2_sym.cptr(), W11aa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W12(c0,i,a2,a) += (    1.00000000) T2(c0,a0,a1,a) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (    1.00000000) V2(c0,w,j,a2) W12(c0,i,a2,a) 
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
    orz::DTensor W12aa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no3_x0_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO3_X0_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, T2b.cptr(), W12aa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_coov_no3_x1_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO3_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, V2_sym.cptr(), W12aa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W14(w,c0,i,a0) += (    1.00000000) V2(w,a2,c0,a1) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (    1.00000000) T2(c0,a0,j,a) W14(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W14caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ocov_coov_no4_x0_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO4_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W14caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no4_x1_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO4_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W14caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W15(w,c0,i,a0) += (    1.00000000) V2(w,c0,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(i,w,j,a) += (    1.00000000) T2(c0,a0,j,a) W15(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W15caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ocov_coov_no5_x0_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO5_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W15caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no5_x1_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO5_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W15caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W17(w,a0,a2,a) += (    1.00000000) V2(w,a2,c0,a1) T2(c0,a0,a1,a) 
  // |-- [    1] --| S2(i,w,j,a) += (   -2.00000000) D2(j,i,a0,a2) W17(w,a0,a2,a) 
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
    orz::DTensor W17aa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no6_x0_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO6_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), W17aa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_coov_no6_x1_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO6_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, W17aa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W18(w,a0,a2,a) += (    1.00000000) V2(w,c0,a2,a1) T2(c0,a0,a1,a) 
  // |-- [    1] --| S2(i,w,j,a) += (    1.00000000) D2(j,i,a0,a2) W18(w,a0,a2,a) 
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
    orz::DTensor W18aa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no7_x0_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO7_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), W18aa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_coov_no7_x1_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO7_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, W18aa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W20(w,a0,j,a) += (    1.00000000) V2(w,j,c0,a1) T2(c0,a0,a1,a) 
  // |-- [    1] --| S2(i,w,j,a) += (   -2.00000000) D1(i,a0) W20(w,a0,j,a) 
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
    orz::DTensor W20aa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no8_x0_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO8_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), W20aa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_coov_no8_x1_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO8_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, W20aa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W22(w,a0,j,a) += (    1.00000000) V2(w,c0,j,a1) T2(c0,a0,a1,a) 
  // |-- [    1] --| S2(i,w,j,a) += (    1.00000000) D1(i,a0) W22(w,a0,j,a) 
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
    orz::DTensor W22aa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no9_x0_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO9_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), W22aa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_coov_no9_x1_type1_eri_c,G_IF_SIGMA_OCOV_COOV_NO9_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, W22aa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // --  Title : sigma_ocov_coov
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W16aa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
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
  // -- Title : sigma_ocov_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W13(i,a1,a0,j) += (    1.00000000) V2(j,a3,a4,a2) D3(i,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) T2(w,a0,a1,a) W13(i,a1,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W13aaa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ocov_coov_no0_x0_type1_eri_o,G_IF_SIGMA_OCOV_COOV_NO0_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W13aaa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no0_x1_type1_eri_o,G_IF_SIGMA_OCOV_COOV_NO0_X1_TYPE1_ERI_O)
      (sa, ia, sj, ij, T2b.cptr(), W13aaa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W16(i,a0) += (    1.00000000) V2(a4,a2,a3,a1) D3(i,a0,a4,a2,a3,a1) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_ocov_coov_no1_x0_type1_eri_o,G_IF_SIGMA_OCOV_COOV_NO1_X0_TYPE1_ERI_O)
    (sa4, ia4, V2_sym.cptr(), W16aa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W19(j,a0,i,a1) += (    1.00000000) V2(a1,a3,a4,a2) D3(j,i,a4,a2,a0,a3) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) T2(w,a0,a1,a) W19(j,a0,i,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W19aaa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa1));
  FC_FUNC(g_if_sigma_ocov_coov_no2_x0_type1_eri_o,G_IF_SIGMA_OCOV_COOV_NO2_X0_TYPE1_ERI_O)
    (sa1, ia1, V2_sym.cptr(), W19aaa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no2_x1_type1_eri_o,G_IF_SIGMA_OCOV_COOV_NO2_X1_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), W19aaa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W21(i,a0,j,a1) += (    1.00000000) V2(j,a2,a3,a1) D2(i,a2,a3,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) T2(w,a0,a1,a) W21(i,a0,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W21aaa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ocov_coov_no3_x0_type1_eri_o,G_IF_SIGMA_OCOV_COOV_NO3_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W21aaa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no3_x1_type1_eri_o,G_IF_SIGMA_OCOV_COOV_NO3_X1_TYPE1_ERI_O)
      (sa, ia, sj, ij, T2b.cptr(), W21aaa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W23(i,a0,j,a2) += (    1.00000000) V2(j,a2,a3,a1) D2(i,a0,a3,a1) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) T2(w,a0,a2,a) W23(i,a0,j,a2) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W23aaa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ocov_coov_no4_x0_type1_eri_o,G_IF_SIGMA_OCOV_COOV_NO4_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W23aaa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no4_x1_type1_eri_o,G_IF_SIGMA_OCOV_COOV_NO4_X1_TYPE1_ERI_O)
      (sa, ia, sj, ij, T2b.cptr(), W23aaa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ocov_coov
  //*-- Entering to take the type 2 contractions --*//

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(i,w,j,a) += (   -0.50000000) T2(w,a0,j,a) W16(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no0_x0_type2_eri_o,G_IF_SIGMA_OCOV_COOV_NO0_X0_TYPE2_ERI_O)
      (sa, ia, T2b.cptr(), W16aa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ocov_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W24(w,a1,a0,a) += (    1.00000000) V2(a,w,v0,c0) T2(c0,a1,a0,v0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -2.00000000) D2(j,i,a1,a0) W24(w,a1,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W24caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ocov_coov_no0_x0_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO0_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W24caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ocov_coov_no0_x1_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO0_X1_TYPE1_ERI_V)
    (sa, ia, W24caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W25(w,a1,a0,a3,a2,a) += (    1.00000000) V2(a,a3,v0,a2) T2(w,a1,a0,v0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) D3(j,i,a3,a0,a1,a2) W25(w,a1,a0,a3,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W25caaaa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ocov_coov_no1_x0_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO1_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W25caaaa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ocov_coov_no1_x1_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO1_X1_TYPE1_ERI_V)
    (sa, ia, W25caaaa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W26(w,a1,a0,a) += (    1.00000000) V2(a,v0,w,c0) T2(c0,a1,a0,v0) 
  // |-- [    1] --| S2(i,w,j,a) += (    1.00000000) D2(j,i,a1,a0) W26(w,a1,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W26caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ocov_coov_no2_x0_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO2_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W26caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ocov_coov_no2_x1_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO2_X1_TYPE1_ERI_V)
    (sa, ia, W26caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W27(w,a1,a0,a3,a2,a) += (    1.00000000) V2(a,v0,a3,a2) T2(w,a1,a0,v0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) D3(j,i,a3,a2,a1,a0) W27(w,a1,a0,a3,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W27caaaa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ocov_coov_no3_x0_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO3_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W27caaaa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ocov_coov_no3_x1_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO3_X1_TYPE1_ERI_V)
    (sa, ia, W27caaaa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W28(w,i,a2,v0) += (    1.00000000) T2(w,a0,a1,v0) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) V2(v0,a,j,a2) W28(w,i,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W28caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ocov_coov_no4_x0_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO4_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W28caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no4_x1_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO4_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W28caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W29(w,a0,j,a) += (    1.00000000) V2(a,w,v0,c0) T2(c0,a0,j,v0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -2.00000000) D1(i,a0) W29(w,a0,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W29caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ocov_coov_no5_x0_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO5_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W29caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ocov_coov_no5_x1_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO5_X1_TYPE1_ERI_V)
    (sa, ia, W29caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W30(i,a0,v0,a) += (    1.00000000) V2(v0,a1,a2,a) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) T2(w,a0,j,v0) W30(i,a0,v0,a) 
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
    orz::DTensor W30aa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ocov_coov_no6_x0_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO6_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W30aa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_coov_no6_x1_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO6_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W30aa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W31(w,a0,j,a) += (    1.00000000) V2(a,v0,w,c0) T2(c0,a0,j,v0) 
  // |-- [    1] --| S2(i,w,j,a) += (    1.00000000) D1(i,a0) W31(w,a0,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W31caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ocov_coov_no7_x0_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO7_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W31caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ocov_coov_no7_x1_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO7_X1_TYPE1_ERI_V)
    (sa, ia, W31caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W32(i,a0,v0,a) += (    1.00000000) V2(v0,a,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) T2(w,a0,j,v0) W32(i,a0,v0,a) 
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
    orz::DTensor W32aa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ocov_coov_no8_x0_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO8_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W32aa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_coov_no8_x1_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO8_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W32aa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W33(w,i,a2,v0) += (    1.00000000) T2(w,a0,a1,v0) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) V2(v0,j,a2,a) W33(w,i,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W33caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ocov_coov_no9_x0_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO9_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W33caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_coov_no9_x1_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO9_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W33caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W34(w,a0,a1,a) += (    1.00000000) V2(a,a2,v0,a1) T2(w,a0,a2,v0) 
  // |-- [    1] --| S2(i,w,j,a) += (    2.00000000) D2(j,i,a0,a1) W34(w,a0,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W34caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ocov_coov_no10_x0_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO10_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W34caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ocov_coov_no10_x1_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO10_X1_TYPE1_ERI_V)
    (sa, ia, W34caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W35(w,a0,a2,a) += (    1.00000000) V2(a,v0,a2,a1) T2(w,a0,a1,v0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) D2(j,i,a0,a2) W35(w,a0,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W35caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ocov_coov_no11_x0_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO11_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W35caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ocov_coov_no11_x1_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO11_X1_TYPE1_ERI_V)
    (sa, ia, W35caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W36(w,a0,j,a) += (    1.00000000) V2(a,v0,j,a1) T2(w,a0,a1,v0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) D1(i,a0) W36(w,a0,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W36caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ocov_coov_no12_x0_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO12_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W36caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ocov_coov_no12_x1_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO12_X1_TYPE1_ERI_V)
    (sa, ia, W36caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W37(w,a0,j,a) += (    1.00000000) V2(a,a1,v0,j) T2(w,a0,a1,v0) 
  // |-- [    1] --| S2(i,w,j,a) += (    2.00000000) D1(i,a0) W37(w,a0,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W37caa_sigma_ocov_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ocov_coov_no13_x0_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO13_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W37caa_sigma_ocov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ocov_coov_no13_x1_type1_eri_v,G_IF_SIGMA_OCOV_COOV_NO13_X1_TYPE1_ERI_V)
    (sa, ia, W37caa_sigma_ocov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ocov_coov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
