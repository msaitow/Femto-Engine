                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ccoo_ccov.h>                                  
                                                                                
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
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccoo_ccov(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(x,w,a1,a0) += (    1.00000000) T2(x,w,a0,v0) Fc1(v0,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) D2(j,a0,i,a1) W0(x,w,a1,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W0ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccoo_ccov_no0_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO0_X0_TYPE0_NOERI)
      (sv0, iv0, T2b.cptr(), W0ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no0_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO0_X1_TYPE0_NOERI)
      (sj, ij, W0ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(x,w,a1,a0) += (    1.00000000) T2(x,w,v0,a0) Fc1(v0,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) D2(j,a1,i,a0) W1(x,w,a1,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W1cca_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no1_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO1_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W1cca_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccov_no1_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO1_X1_TYPE0_NOERI)
        (sa0, ia0, sj, ij, W1cca_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(j,v0) += (    1.00000000) D1(j,a0) Fc1(v0,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,w,v0,i) W2(j,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W2v_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xv(symblockinfo, sj));
    FC_FUNC(g_if_sigma_ccoo_ccov_no2_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO2_X0_TYPE0_NOERI)
      (sj, ij, W2v_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    for(int si = 0;si < nir;++si){ 
    for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
      T2b = T2.get_amp2(ii);
      FC_FUNC(g_if_sigma_ccoo_ccov_no2_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO2_X1_TYPE0_NOERI)
        (si, ii, sj, ij, T2b.cptr(), W2v_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    } // End ii
    } // End si
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(j,v0) += (    1.00000000) D1(j,a0) Fc1(v0,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,i,v0) W3(j,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      double W3_sigma_ccoo_ccov(0);
      FC_FUNC(g_if_sigma_ccoo_ccov_no3_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO3_X0_TYPE0_NOERI)
        (sj, ij, sv0, iv0, &W3_sigma_ccoo_ccov, nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccov_no3_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO3_X1_TYPE0_NOERI)
        (sj, ij, sv0, iv0, T2b.cptr(), &W3_sigma_ccoo_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(x,w,j,v0) += (    1.00000000) T2(x,w,a0,v0) D1(j,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) Fc1(v0,i) W4(x,w,j,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      orz::DTensor W4cc_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sj^sv0));
      FC_FUNC(g_if_sigma_ccoo_ccov_no4_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO4_X0_TYPE0_NOERI)
        (sj, ij, sv0, iv0, T2b.cptr(), W4cc_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
      FC_FUNC(g_if_sigma_ccoo_ccov_no4_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO4_X1_TYPE0_NOERI)
        (sj, ij, sv0, iv0, W4cc_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W5(x,w,j,v0) += (    1.00000000) T2(x,w,v0,a0) D1(j,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) Fc1(v0,i) W5(x,w,j,v0) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W5ccv_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccv(symblockinfo, sj));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_ccov_no5_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO5_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W5ccv_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no5_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO5_X1_TYPE0_NOERI)
      (sj, ij, W5ccv_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W6(i,v0) += (    1.00000000) D1(i,a0) Fc1(v0,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,v0,j) W6(i,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W6av_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_ccov_no6_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO6_X0_TYPE0_NOERI)
    (W6av_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no6_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO6_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W6av_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W7(i,v0) += (    1.00000000) D1(i,a0) Fc1(v0,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,x,v0,j) W7(i,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W7av_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_ccov_no7_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO7_X0_TYPE0_NOERI)
    (W7av_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no7_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO7_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W7av_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W8(x,w,i,v0) += (    1.00000000) T2(x,w,a0,v0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) Fc1(v0,j) W8(x,w,i,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    orz::DTensor W8cca_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no8_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO8_X0_TYPE0_NOERI)
      (sv0, iv0, T2b.cptr(), W8cca_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccov_no8_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO8_X1_TYPE0_NOERI)
        (sj, ij, sv0, iv0, W8cca_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W9(x,w,i,v0) += (    1.00000000) T2(x,w,v0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) Fc1(v0,j) W9(x,w,i,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9ccav_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccav(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_ccov_no9_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO9_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W9ccav_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no9_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO9_X1_TYPE0_NOERI)
      (sj, ij, W9ccav_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) T2(x,w,v0,i) Fc1(v0,j) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccov_no10_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO10_X0_TYPE0_NOERI)
        (si, ii, sj, ij, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) T2(x,w,i,v0) Fc1(v0,j) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccov_no11_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO11_X0_TYPE0_NOERI)
        (sj, ij, sv0, iv0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) T2(w,x,v0,j) Fc1(v0,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no12_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO12_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) T2(x,w,v0,j) Fc1(v0,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no13_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOV_NO13_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
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
  // -- Title : sigma_ccoo_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W48(w,c0,i,v0) += (    1.00000000) V2(w,a0,v0,c0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(c0,x,v0,j) W48(w,c0,i,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W48cav_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccov_no0_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO0_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W48cav_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no0_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO0_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W48cav_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W49(x,c0,i,v0) += (    1.00000000) V2(x,a0,v0,c0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(c0,w,v0,j) W49(x,c0,i,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W49cav_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccov_no1_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO1_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W49cav_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no1_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO1_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W49cav_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W50(w,c0,i,v0) += (    1.00000000) V2(w,a0,v0,c0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,c0,v0,j) W50(w,c0,i,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W50cav_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccov_no2_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO2_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W50cav_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no2_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO2_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W50cav_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W51(x,c0,i,v0) += (    1.00000000) V2(x,a0,v0,c0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) T2(w,c0,v0,j) W51(x,c0,i,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W51cav_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccov_no3_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO3_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W51cav_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no3_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO3_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W51cav_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W52(w,c0,i,v0) += (    1.00000000) V2(w,c0,v0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(c0,x,v0,j) W52(w,c0,i,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W52cav_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccov_no4_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO4_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W52cav_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no4_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO4_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W52cav_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W53(x,c0,i,v0) += (    1.00000000) V2(x,c0,v0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) T2(c0,w,v0,j) W53(x,c0,i,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W53cav_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccov_no5_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO5_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W53cav_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no5_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO5_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W53cav_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W54(w,c0,i,v0) += (    1.00000000) V2(w,c0,v0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) T2(x,c0,v0,j) W54(w,c0,i,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W54cav_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccov_no6_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO6_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W54cav_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no6_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO6_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W54cav_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W55(x,c0,i,v0) += (    1.00000000) V2(x,c0,v0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(w,c0,v0,j) W55(x,c0,i,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W55cav_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccov_no7_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO7_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W55cav_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no7_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO7_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W55cav_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    8.00000000) V2(w,i,v0,c0) T2(c0,x,v0,j) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no8_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO8_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) V2(x,i,v0,c0) T2(c0,w,v0,j) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no9_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO9_X0_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) V2(w,i,v0,c0) T2(x,c0,v0,j) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no10_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO10_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) V2(x,i,v0,c0) T2(w,c0,v0,j) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no11_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO11_X0_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) V2(x,c0,v0,i) T2(w,c0,v0,j) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no12_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO12_X0_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) V2(w,c0,v0,i) T2(x,c0,v0,j) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no13_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO13_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) V2(x,c0,v0,i) T2(c0,w,v0,j) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no14_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO14_X0_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) V2(w,c0,v0,i) T2(c0,x,v0,j) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no15_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOV_NO15_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ccoo_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W42(w,x,a0,j) += (    1.00000000) V2(j,x,v0,c0) T2(w,c0,v0,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) D1(i,a0) W42(w,x,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W42cc_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_ccov_no0_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO0_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W42cc_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no0_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO0_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, W42cc_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W43(x,w,a0,j) += (    1.00000000) V2(j,w,v0,c0) T2(x,c0,v0,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) D1(i,a0) W43(x,w,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W43cc_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_ccov_no1_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO1_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W43cc_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no1_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO1_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, W43cc_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W44(x,w,a0,j) += (    1.00000000) V2(j,w,v0,c0) T2(x,c0,a0,v0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) D1(i,a0) W44(x,w,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W44cca_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccoo_ccov_no2_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO2_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W44cca_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ccoo_ccov_no2_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO2_X1_TYPE1_ERI_O)
    (sj, ij, W44cca_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W45(w,x,a0,j) += (    1.00000000) V2(j,x,v0,c0) T2(w,c0,a0,v0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) D1(i,a0) W45(w,x,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W45cca_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccoo_ccov_no3_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO3_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W45cca_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ccoo_ccov_no3_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO3_X1_TYPE1_ERI_O)
    (sj, ij, W45cca_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W58(x,w,a0,j) += (    1.00000000) V2(j,v0,w,c0) T2(x,c0,v0,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) D1(i,a0) W58(x,w,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W58cc_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_ccov_no4_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO4_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W58cc_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no4_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO4_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, W58cc_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W59(w,x,a0,j) += (    1.00000000) V2(j,v0,x,c0) T2(w,c0,v0,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) D1(i,a0) W59(w,x,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W59cc_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_ccov_no5_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO5_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W59cc_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no5_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO5_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, W59cc_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W60(x,w,a0,j) += (    1.00000000) V2(j,v0,w,c0) T2(x,c0,a0,v0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) D1(i,a0) W60(x,w,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W60cca_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccoo_ccov_no6_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO6_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W60cca_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ccoo_ccov_no6_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO6_X1_TYPE1_ERI_O)
    (sj, ij, W60cca_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W61(w,x,a0,j) += (    1.00000000) V2(j,v0,x,c0) T2(w,c0,a0,v0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) D1(i,a0) W61(w,x,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W61cca_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccoo_ccov_no7_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO7_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W61cca_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ccoo_ccov_no7_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO7_X1_TYPE1_ERI_O)
    (sj, ij, W61cca_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    8.00000000) V2(j,x,v0,c0) T2(w,c0,i,v0) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccoo_ccov_no8_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO8_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) V2(j,w,v0,c0) T2(x,c0,i,v0) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccoo_ccov_no9_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO9_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) V2(j,w,v0,c0) T2(x,c0,v0,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_ccov_no10_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO10_X0_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) V2(j,x,v0,c0) T2(w,c0,v0,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_ccov_no11_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOV_NO11_X0_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@type(2).declaration(begin)
  // --  Title : sigma_ccoo_ccov
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W10ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W11ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W12ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W13ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W14ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W15ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W16ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W17ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W20ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W21ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W22ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W23ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W36ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W37ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W38ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W39ccaa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
//-@type(2).declaration(end)

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
  // -- Title : sigma_ccoo_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W10(w,x,a0,a1) += (    1.00000000) V2(v0,c0,x,a1) T2(c0,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no0_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO0_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W10ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W11(x,w,a0,a1) += (    1.00000000) V2(v0,c0,w,a1) T2(c0,x,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no1_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO1_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W11ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W12(x,w,a0,a1) += (    1.00000000) V2(v0,c0,w,a1) T2(x,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no2_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO2_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W12ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W13(w,x,a0,a1) += (    1.00000000) V2(v0,c0,x,a1) T2(w,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no3_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO3_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W13ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W14(w,x,a0,a1) += (    1.00000000) V2(v0,a1,x,c0) T2(c0,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no4_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO4_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W14ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W15(x,w,a0,a1) += (    1.00000000) V2(v0,a1,w,c0) T2(c0,x,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no5_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO5_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W15ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W16(x,w,a0,a1) += (    1.00000000) V2(v0,a1,w,c0) T2(x,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no6_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO6_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W16ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W17(w,x,a0,a1) += (    1.00000000) V2(v0,a1,x,c0) T2(w,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no7_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO7_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W17ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W18(j,i,a0,v0) += (    1.00000000) V2(v0,a2,a3,a1) D3(j,a2,i,a0,a3,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,x,a0,v0) W18(j,i,a0,v0) 
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
    orz::DTensor W18aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no8_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO8_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W18aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no8_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO8_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W18aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W19(j,i,a0,v0) += (    1.00000000) V2(v0,a2,a3,a1) D3(j,a0,i,a2,a3,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,a0,v0) W19(j,i,a0,v0) 
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
    orz::DTensor W19aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no9_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO9_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W19aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no9_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO9_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W19aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W20(w,x,a0,i) += (    1.00000000) V2(v0,c0,x,i) T2(c0,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no10_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO10_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W20ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W21(x,w,a0,i) += (    1.00000000) V2(v0,c0,w,i) T2(c0,x,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no11_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO11_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W21ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W22(x,w,a0,i) += (    1.00000000) V2(v0,c0,w,i) T2(x,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no12_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO12_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W22ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W23(w,x,a0,i) += (    1.00000000) V2(v0,c0,x,i) T2(w,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no13_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO13_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W23ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W24(j,a0,i,v0) += (    1.00000000) V2(v0,a2,i,a1) D2(j,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,x,a0,v0) W24(j,a0,i,v0) 
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
    orz::DTensor W24aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no14_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO14_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W24aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no14_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO14_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W24aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W25(j,a0,i,v0) += (    1.00000000) V2(v0,a2,i,a1) D2(j,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,a0,v0) W25(j,a0,i,v0) 
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
    orz::DTensor W25aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no15_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO15_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W25aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no15_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO15_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W25aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W26(w,c0,j,v0) += (    1.00000000) V2(v0,c0,w,a0) D1(j,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,c0,i,v0) W26(w,c0,j,v0) 
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
    orz::DTensor W26cc_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no16_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO16_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W26cc_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no16_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO16_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W26cc_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W27(x,c0,j,v0) += (    1.00000000) V2(v0,c0,x,a0) D1(j,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(w,c0,i,v0) W27(x,c0,j,v0) 
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
    orz::DTensor W27cc_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no17_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO17_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W27cc_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no17_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO17_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W27cc_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W28(w,c0,j,v0) += (    1.00000000) V2(v0,c0,w,a0) D1(j,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) T2(c0,x,i,v0) W28(w,c0,j,v0) 
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
    orz::DTensor W28cc_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no18_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO18_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W28cc_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no18_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO18_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W28cc_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W29(x,c0,j,v0) += (    1.00000000) V2(v0,c0,x,a0) D1(j,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(c0,w,i,v0) W29(x,c0,j,v0) 
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
    orz::DTensor W29cc_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no19_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO19_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W29cc_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no19_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO19_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W29cc_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W30(w,c0,j,v0) += (    1.00000000) V2(v0,a0,w,c0) D1(j,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) T2(x,c0,i,v0) W30(w,c0,j,v0) 
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
    orz::DTensor W30cc_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no20_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO20_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W30cc_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no20_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO20_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W30cc_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W31(x,c0,j,v0) += (    1.00000000) V2(v0,a0,x,c0) D1(j,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(w,c0,i,v0) W31(x,c0,j,v0) 
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
    orz::DTensor W31cc_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no21_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO21_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W31cc_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no21_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO21_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W31cc_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W32(w,c0,j,v0) += (    1.00000000) V2(v0,a0,w,c0) D1(j,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(c0,x,i,v0) W32(w,c0,j,v0) 
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
    orz::DTensor W32cc_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no22_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO22_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W32cc_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no22_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO22_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W32cc_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W33(x,c0,j,v0) += (    1.00000000) V2(v0,a0,x,c0) D1(j,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) T2(c0,w,i,v0) W33(x,c0,j,v0) 
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
    orz::DTensor W33cc_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no23_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO23_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W33cc_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no23_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO23_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W33cc_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W34(j,v0) += (    1.00000000) V2(v0,a1,a2,a0) D2(j,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,i,v0) W34(j,v0) 
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
    double W34_sigma_ccoo_ccov(0);
    FC_FUNC(g_if_sigma_ccoo_ccov_no24_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO24_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), &W34_sigma_ccoo_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no24_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO24_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), &W34_sigma_ccoo_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W35(j,v0) += (    1.00000000) V2(v0,a1,a2,a0) D2(j,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,x,i,v0) W35(j,v0) 
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
    double W35_sigma_ccoo_ccov(0);
    FC_FUNC(g_if_sigma_ccoo_ccov_no25_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO25_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), &W35_sigma_ccoo_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no25_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO25_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), &W35_sigma_ccoo_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W36(x,w,a0,i) += (    1.00000000) V2(v0,i,w,c0) T2(c0,x,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no26_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO26_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W36ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W37(w,x,a0,i) += (    1.00000000) V2(v0,i,x,c0) T2(c0,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no27_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO27_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W37ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W38(x,w,a0,i) += (    1.00000000) V2(v0,i,w,c0) T2(x,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no28_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO28_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W38ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W39(w,x,a0,i) += (    1.00000000) V2(v0,i,x,c0) T2(w,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_ccoo_ccov_no29_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO29_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W39ccaa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W40(j,a0,i,v0) += (    1.00000000) V2(v0,i,a2,a1) D2(j,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,x,a0,v0) W40(j,a0,i,v0) 
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
    orz::DTensor W40aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no30_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO30_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W40aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no30_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO30_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W40aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W41(j,a0,i,v0) += (    1.00000000) V2(v0,i,a2,a1) D2(j,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,w,a0,v0) W41(j,a0,i,v0) 
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
    orz::DTensor W41aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no31_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO31_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W41aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no31_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO31_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W41aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W46(i,a0,j,v0) += (    1.00000000) V2(v0,a2,j,a1) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,x,a0,v0) W46(i,a0,j,v0) 
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
    orz::DTensor W46aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no32_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO32_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W46aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no32_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO32_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W46aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W47(i,a0,j,v0) += (    1.00000000) V2(v0,a2,j,a1) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,a0,v0) W47(i,a0,j,v0) 
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
    orz::DTensor W47aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no33_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO33_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W47aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no33_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO33_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W47aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W56(i,v0) += (    1.00000000) V2(v0,a1,a2,a0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,w,j,v0) W56(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W56a_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccoo_ccov_no34_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO34_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W56a_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no34_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO34_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W56a_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W57(i,v0) += (    1.00000000) V2(v0,a1,a2,a0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,v0,j) W57(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W57a_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccoo_ccov_no35_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO35_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W57a_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no35_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO35_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W57a_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W62(i,a0,j,v0) += (    1.00000000) V2(v0,j,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,x,a0,v0) W62(i,a0,j,v0) 
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
    orz::DTensor W62aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no36_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO36_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W62aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no36_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO36_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W62aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W63(i,a0,j,v0) += (    1.00000000) V2(v0,j,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,a0,v0) W63(i,a0,j,v0) 
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
    orz::DTensor W63aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no37_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO37_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W63aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no37_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO37_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W63aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W64(i,v0) += (    1.00000000) V2(v0,a1,i,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,w,j,v0) W64(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W64a_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccoo_ccov_no38_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO38_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W64a_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no38_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO38_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W64a_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W65(i,v0) += (    1.00000000) V2(v0,a1,i,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,v0,j) W65(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W65a_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccoo_ccov_no39_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO39_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W65a_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no39_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO39_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W65a_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| W66(a0,i,j,v0) += (    1.00000000) V2(v0,j,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,x,a0,v0) W66(a0,i,j,v0) 
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
    orz::DTensor W66aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no40_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO40_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W66aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no40_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO40_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W66aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| W67(a0,i,j,v0) += (    1.00000000) V2(v0,j,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,a0,v0) W67(a0,i,j,v0) 
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
    orz::DTensor W67aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no41_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO41_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W67aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no41_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO41_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W67aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| W68(j,v0) += (    1.00000000) V2(v0,a1,j,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,i,v0) W68(j,v0) 
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
    double W68_sigma_ccoo_ccov(0);
    FC_FUNC(g_if_sigma_ccoo_ccov_no42_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO42_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), &W68_sigma_ccoo_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no42_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO42_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), &W68_sigma_ccoo_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| W69(j,v0) += (    1.00000000) V2(v0,a1,j,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,x,i,v0) W69(j,v0) 
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
    double W69_sigma_ccoo_ccov(0);
    FC_FUNC(g_if_sigma_ccoo_ccov_no43_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO43_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), &W69_sigma_ccoo_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no43_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO43_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), &W69_sigma_ccoo_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) V2(v0,j,x,c0) T2(w,c0,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no44_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO44_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) V2(v0,j,w,c0) T2(x,c0,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no45_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO45_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) V2(v0,j,x,c0) T2(c0,w,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no46_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO46_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) V2(v0,j,w,c0) T2(c0,x,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no47_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO47_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   48] -- 
  // |-- [    0] --| W70(j,v0) += (    1.00000000) V2(v0,j,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(w,x,i,v0) W70(j,v0) 
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
    double W70_sigma_ccoo_ccov(0);
    FC_FUNC(g_if_sigma_ccoo_ccov_no48_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO48_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), &W70_sigma_ccoo_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no48_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO48_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), &W70_sigma_ccoo_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   49] -- 
  // |-- [    0] --| W71(j,v0) += (    1.00000000) V2(v0,j,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,w,i,v0) W71(j,v0) 
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
    double W71_sigma_ccoo_ccov(0);
    FC_FUNC(g_if_sigma_ccoo_ccov_no49_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO49_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), &W71_sigma_ccoo_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no49_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO49_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), &W71_sigma_ccoo_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   50] -- 
  // |-- [    0] --| W72(a0,j,i,v0) += (    1.00000000) V2(v0,i,j,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,x,a0,v0) W72(a0,j,i,v0) 
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
    orz::DTensor W72aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no50_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO50_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W72aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no50_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO50_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W72aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   51] -- 
  // |-- [    0] --| W73(a0,j,i,v0) += (    1.00000000) V2(v0,i,j,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,w,a0,v0) W73(a0,j,i,v0) 
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
    orz::DTensor W73aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no51_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO51_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W73aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no51_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO51_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W73aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   52] -- 
  // |-- [    0] --| W74(i,v0) += (    1.00000000) V2(v0,i,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(x,w,j,v0) W74(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W74a_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccoo_ccov_no52_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO52_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W74a_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no52_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO52_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W74a_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   53] -- 
  // |-- [    0] --| W75(i,v0) += (    1.00000000) V2(v0,i,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,w,v0,j) W75(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W75a_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccoo_ccov_no53_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO53_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W75a_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no53_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO53_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W75a_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   54] -- 
  // |-- [    0] --| W76(j,i,a0,v0) += (    1.00000000) V2(v0,a1,a2,a0) D2(j,a2,i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,a0,v0) W76(j,i,a0,v0) 
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
    orz::DTensor W76aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no54_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO54_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W76aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no54_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO54_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W76aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   55] -- 
  // |-- [    0] --| W77(j,i,a0,v0) += (    1.00000000) V2(v0,a1,a2,a0) D2(j,a1,i,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,x,a0,v0) W77(j,i,a0,v0) 
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
    orz::DTensor W77aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no55_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO55_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W77aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no55_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO55_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W77aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   56] -- 
  // |-- [    0] --| W78(j,i,a0,v0) += (    1.00000000) V2(v0,a1,i,a0) D1(j,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,a0,v0) W78(j,i,a0,v0) 
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
    orz::DTensor W78aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no56_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO56_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W78aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no56_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO56_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W78aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   57] -- 
  // |-- [    0] --| W79(j,i,a0,v0) += (    1.00000000) V2(v0,a1,i,a0) D1(j,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,x,a0,v0) W79(j,i,a0,v0) 
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
    orz::DTensor W79aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no57_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO57_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W79aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no57_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO57_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W79aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   58] -- 
  // |-- [    0] --| W80(j,i,a0,v0) += (    1.00000000) V2(v0,i,a1,a0) D1(j,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,x,a0,v0) W80(j,i,a0,v0) 
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
    orz::DTensor W80aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no58_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO58_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W80aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no58_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO58_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W80aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   59] -- 
  // |-- [    0] --| W81(j,i,a0,v0) += (    1.00000000) V2(v0,i,a1,a0) D1(j,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,w,a0,v0) W81(j,i,a0,v0) 
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
    orz::DTensor W81aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no59_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO59_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W81aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no59_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO59_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W81aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   60] -- 
  // |-- [    0] --| W82(i,j,a0,v0) += (    1.00000000) V2(v0,a1,j,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,w,a0,v0) W82(i,j,a0,v0) 
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
    orz::DTensor W82aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no60_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO60_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W82aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no60_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO60_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W82aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   61] -- 
  // |-- [    0] --| W83(i,j,a0,v0) += (    1.00000000) V2(v0,a1,j,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,x,a0,v0) W83(i,j,a0,v0) 
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
    orz::DTensor W83aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no61_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO61_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W83aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no61_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO61_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W83aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   62] -- 
  // |-- [    0] --| W84(i,j,a0,v0) += (    1.00000000) V2(v0,j,a1,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,x,a0,v0) W84(i,j,a0,v0) 
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
    orz::DTensor W84aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no62_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO62_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W84aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no62_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO62_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W84aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   63] -- 
  // |-- [    0] --| W85(i,j,a0,v0) += (    1.00000000) V2(v0,j,a1,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,a0,v0) W85(i,j,a0,v0) 
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
    orz::DTensor W85aa_sigma_ccoo_ccov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_ccoo_ccov_no63_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO63_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W85aa_sigma_ccoo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccov_no63_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO63_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W85aa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   64] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) V2(v0,j,i,a0) T2(w,x,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no64_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO64_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   65] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) V2(v0,j,i,a0) T2(x,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no65_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO65_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   66] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) V2(v0,i,j,a0) T2(x,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no66_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO66_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   67] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) V2(v0,i,j,a0) T2(w,x,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no67_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCOV_NO67_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@type(2).contraction(begin)
  // -- Title : sigma_ccoo_ccov
  //*-- Entering to take the type 2 contractions --*//

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -1.00000000) D2(j,a1,i,a0) W10(w,x,a0,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no0_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO0_X0_TYPE2_ERI_V)
      (sj, ij, W10ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -1.00000000) D2(j,a0,i,a1) W11(x,w,a0,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no1_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO1_X0_TYPE2_ERI_V)
      (sj, ij, W11ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) D2(j,a0,i,a1) W12(x,w,a0,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no2_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO2_X0_TYPE2_ERI_V)
      (sj, ij, W12ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) D2(j,a1,i,a0) W13(w,x,a0,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no3_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO3_X0_TYPE2_ERI_V)
      (sj, ij, W13ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -1.00000000) D2(j,a0,i,a1) W14(w,x,a0,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no4_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO4_X0_TYPE2_ERI_V)
      (sj, ij, W14ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -1.00000000) D2(j,a1,i,a0) W15(x,w,a0,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no5_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO5_X0_TYPE2_ERI_V)
      (sj, ij, W15ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -1.00000000) D2(j,a0,i,a1) W16(x,w,a0,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no6_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO6_X0_TYPE2_ERI_V)
      (sj, ij, W16ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -1.00000000) D2(j,a1,i,a0) W17(w,x,a0,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no7_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO7_X0_TYPE2_ERI_V)
      (sj, ij, W17ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -1.00000000) D1(j,a0) W20(w,x,a0,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no8_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO8_X0_TYPE2_ERI_V)
      (sj, ij, W20ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) D1(j,a0) W21(x,w,a0,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no9_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO9_X0_TYPE2_ERI_V)
      (sj, ij, W21ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) D1(j,a0) W22(x,w,a0,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no10_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO10_X0_TYPE2_ERI_V)
      (sj, ij, W22ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) D1(j,a0) W23(w,x,a0,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no11_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO11_X0_TYPE2_ERI_V)
      (sj, ij, W23ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -1.00000000) D1(j,a0) W36(x,w,a0,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no12_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO12_X0_TYPE2_ERI_V)
      (sj, ij, W36ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) D1(j,a0) W37(w,x,a0,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no13_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO13_X0_TYPE2_ERI_V)
      (sj, ij, W37ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) D1(j,a0) W38(x,w,a0,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no14_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO14_X0_TYPE2_ERI_V)
      (sj, ij, W38ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -1.00000000) D1(j,a0) W39(w,x,a0,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccov_no15_x0_type2_eri_v,G_IF_SIGMA_CCOO_CCOV_NO15_X0_TYPE2_ERI_V)
      (sj, ij, W39ccaa_sigma_ccoo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
//-@type(2).contraction(end)
  orz::world().barrier();

  } // End Femto
  //*-- FEMTO ends --//*

//-@loadERI(v,end)

//-@ERI.contractions(end)

//-@D4C.contractions(begin)

//-@D4C.contractions(end)

  // Do timing!
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccoo_ccov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
