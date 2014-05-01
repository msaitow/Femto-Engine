                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ooov_ooov.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//       #                #########      #     #   # 
//  ########## ##########         #   #######  #   # 
//      #    #         #          #    # #     #   # 
//      #    #        #   ########     # #     #   # 
//     #     #     # #           #  ##########    #  
//    #   # #       #            #       #       #   
//   #     #         #    ########       #     ##    

//                                   Generated date : Sun Apr 20 10:26:05 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ooov_ooov(const orz::mr::Input &ctinp,                                    
                                  const orz::mr::SymBlockInfo &symblockinfo,                                 
                                  const orz::mr::HintMO &hintmo,                                             
                                  const int alloc_type,                                                      
                                  const orz::mr::RdmPack &rdmPack,                                           
                                  const orz::DTensor &rdm4,                                                  
                                  const orz::mr::BareAmpPack &T2,                             
                                  const int num_sigma,
                                  const double E0) {
                                                                                                                 
                                                                                                                 
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
  // |-- [    0] --| W0(k,a3,a2,i,j,a1) += (    1.00000000) D3(k,i,a3,a0,a2,j) Fc1(a1,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) T2(a2,a3,a,a1) W0(k,a3,a2,i,j,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W0aaaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sa1));
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO0_X0_TYPE0_NOERI)
      (sa1, ia1, W0aaaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x1_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO0_X1_TYPE0_NOERI)
        (sa, ia, sa1, ia1, T2b.cptr(), W0aaaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
      retval.acc_amp2(ia, S2b);
    } // End ia
    } // End sa
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(j,i,a0,a) += (    1.00000000) T2(a2,a1,a0,a) D2(j,a1,i,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) Fc1(k,a0) W1(j,i,a0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W1aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO1_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W1aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x1_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO1_X1_TYPE0_NOERI)
      (sa, ia, W1aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(k,i,j,v0) += (    1.00000000) T2(a2,a1,v0,a0) D3(k,i,a2,j,a1,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) Fc1(v0,a) W2(k,i,j,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W2aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no2_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO2_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W2aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no2_x1_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO2_X1_TYPE0_NOERI)
      (sa, ia, W2aaav_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(j,i,k,v0) += (    1.00000000) T2(a1,a0,v0,k) D2(j,a1,i,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) Fc1(v0,a) W3(j,i,k,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sk = 0;sk < nir;++sk){ 
  for(int ik = symblockinfo.psym()(sk,I_O,I_BEGIN);ik <= symblockinfo.psym()(sk,I_O,I_END);++ik){ 
    T2b = T2.get_amp2(ik);
    orz::DTensor W3aav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sk));
    FC_FUNC(g_if_sigma_ooov_ooov_no3_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO3_X0_TYPE0_NOERI)
      (sk, ik, T2b.cptr(), W3aav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_ooov_ooov_no3_x1_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO3_X1_TYPE0_NOERI)
        (sa, ia, sk, ik, W3aav_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
      retval.acc_amp2(ia, S2b);
    } // End ia
    } // End sa
  // --> @[k, "active"] [notNeeded]
  } // End ik
  } // End sk
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(k,a2,i,j,a1,a3) += (    1.00000000) D3(k,i,a2,j,a0,a1) Fc1(a3,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) T2(a2,a3,a,a1) W4(k,a2,i,j,a1,a3) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W4aaaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sa1));
    FC_FUNC(g_if_sigma_ooov_ooov_no4_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO4_X0_TYPE0_NOERI)
      (sa1, ia1, W4aaaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_ooov_ooov_no4_x1_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO4_X1_TYPE0_NOERI)
        (sa, ia, sa1, ia1, T2b.cptr(), W4aaaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
      retval.acc_amp2(ia, S2b);
    } // End ia
    } // End sa
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W5(j,i,a1,a2) += (    1.00000000) D2(j,a1,i,a0) Fc1(a2,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) T2(a2,a1,k,a) W5(j,i,a1,a2) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W5aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ooov_ooov_no5_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO5_X0_TYPE0_NOERI)
    (W5aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no5_x1_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO5_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W5aaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W6(k,a3,i,a1,j,a2) += (    1.00000000) D3(k,i,a3,a1,a0,j) Fc1(a2,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) T2(a2,a3,a,a1) W6(k,a3,i,a1,j,a2) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W6aaaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sa1));
    FC_FUNC(g_if_sigma_ooov_ooov_no6_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO6_X0_TYPE0_NOERI)
      (sa1, ia1, W6aaaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_ooov_ooov_no6_x1_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO6_X1_TYPE0_NOERI)
        (sa, ia, sa1, ia1, T2b.cptr(), W6aaaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
      retval.acc_amp2(ia, S2b);
    } // End ia
    } // End sa
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W7(j,i,a2,a1) += (    1.00000000) D2(j,a0,i,a2) Fc1(a1,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) T2(a2,a1,k,a) W7(j,i,a2,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W7aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ooov_ooov_no7_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO7_X0_TYPE0_NOERI)
    (W7aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no7_x1_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO7_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W7aaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(i,j,k,a) += (    1.00000000) E0 T2(a1,a2,a0,a) D3(k,i,a2,j,a1,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no8_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO8_X0_TYPE0_NOERI)
      (sa, ia, &E0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(i,j,k,a) += (    1.00000000) E0 T2(a0,a1,k,a) D2(j,a1,i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no9_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO9_X0_TYPE0_NOERI)
      (sa, ia, &E0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W304(j,i,a0,a) += (    1.00000000) T2(a3,a2,a1,a) D3(j,a2,i,a0,a1,a3) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) Fc1(k,a0) W304(j,i,a0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W304aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no10_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO10_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W304aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no10_x1_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO10_X1_TYPE0_NOERI)
      (sa, ia, W304aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W305(k,a0,j,a) += (    1.00000000) T2(a3,a2,a1,a) D3(k,a0,a3,a1,a2,j) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) Fc1(i,a0) W305(k,a0,j,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W305aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no11_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO11_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W305aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no11_x1_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO11_X1_TYPE0_NOERI)
      (sa, ia, W305aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W306(j,a0,k,a) += (    1.00000000) T2(a2,a1,k,a) D2(j,a1,a0,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) Fc1(i,a0) W306(j,a0,k,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W306aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no12_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO12_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W306aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no12_x1_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO12_X1_TYPE0_NOERI)
      (sa, ia, W306aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W307(k,i,a0,a) += (    1.00000000) T2(a3,a2,a1,a) D3(k,i,a3,a1,a2,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) Fc1(j,a0) W307(k,i,a0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W307aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no13_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO13_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W307aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no13_x1_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO13_X1_TYPE0_NOERI)
      (sa, ia, W307aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W308(i,a0,k,a) += (    1.00000000) T2(a2,a1,k,a) D2(i,a2,a0,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) Fc1(j,a0) W308(i,a0,k,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W308aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no14_x0_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO14_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W308aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no14_x1_type0_noeri,G_IF_SIGMA_OOOV_OOOV_NO14_X1_TYPE0_NOERI)
      (sa, ia, W308aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  } // End Femto
  //*-- FEMTO ends --//*

//-@ERI.contractions(begin)

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
  // -- Title : sigma_ooov_ooov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W8(j,i,a1,a0,k,a2) += (    1.00000000) V2(k,a3,a4,a2) D3(j,a1,i,a3,a4,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) T2(a0,a1,a2,a) W8(j,i,a1,a0,k,a2) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W8aaaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sk));
  FC_FUNC(g_if_sigma_ooov_ooov_no0_x0_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO0_X0_TYPE1_ERI_O)
    (sk, ik, V2_sym.cptr(), W8aaaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x1_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO0_X1_TYPE1_ERI_O)
      (sa, ia, sk, ik, T2b.cptr(), W8aaaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W9(j,i,a1,a0,k,a3) += (    1.00000000) V2(k,a3,a4,a2) D3(j,a1,i,a0,a4,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) T2(a0,a1,a3,a) W9(j,i,a1,a0,k,a3) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9aaaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sk));
  FC_FUNC(g_if_sigma_ooov_ooov_no1_x0_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO1_X0_TYPE1_ERI_O)
    (sk, ik, V2_sym.cptr(), W9aaaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x1_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO1_X1_TYPE1_ERI_O)
      (sa, ia, sk, ik, T2b.cptr(), W9aaaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W298(j,i,a1,a0,k,a2) += (    1.00000000) V2(k,a3,a4,a2) D3(j,a0,i,a3,a1,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) T2(a2,a0,a1,a) W298(j,i,a1,a0,k,a2) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W298aaaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sk));
  FC_FUNC(g_if_sigma_ooov_ooov_no2_x0_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO2_X0_TYPE1_ERI_O)
    (sk, ik, V2_sym.cptr(), W298aaaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no2_x1_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO2_X1_TYPE1_ERI_O)
      (sa, ia, sk, ik, T2b.cptr(), W298aaaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W299(j,i,a0,a1) += (    1.00000000) V2(a1,a3,a4,a2) D3(j,a0,i,a3,a4,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) T2(a1,a0,k,a) W299(j,i,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W299aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa1));
  FC_FUNC(g_if_sigma_ooov_ooov_no3_x0_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO3_X0_TYPE1_ERI_O)
    (sa1, ia1, V2_sym.cptr(), W299aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no3_x1_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO3_X1_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), W299aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W300(j,i,a0,a1) += (    1.00000000) V2(a1,a3,a4,a2) D3(j,a3,i,a0,a4,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) T2(a0,a1,k,a) W300(j,i,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W300aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa1));
  FC_FUNC(g_if_sigma_ooov_ooov_no4_x0_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO4_X0_TYPE1_ERI_O)
    (sa1, ia1, V2_sym.cptr(), W300aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no4_x1_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO4_X1_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), W300aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W301(j,i,a1,a0,k,a2) += (    1.00000000) V2(k,a3,a4,a2) D3(j,a4,i,a3,a1,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) T2(a0,a2,a1,a) W301(j,i,a1,a0,k,a2) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W301aaaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sk));
  FC_FUNC(g_if_sigma_ooov_ooov_no5_x0_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO5_X0_TYPE1_ERI_O)
    (sk, ik, V2_sym.cptr(), W301aaaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no5_x1_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO5_X1_TYPE1_ERI_O)
      (sa, ia, sk, ik, T2b.cptr(), W301aaaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W302(a0,a4,a3,a) += (    1.00000000) V2(a4,a2,a3,a1) T2(a1,a2,a0,a) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) D3(k,i,a4,j,a3,a0) W302(a0,a4,a3,a) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W302aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa4^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no6_x0_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO6_X0_TYPE1_ERI_O)
      (sa, ia, sa4, ia4, T2b.cptr(), V2_sym.cptr(), W302aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no6_x1_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO6_X1_TYPE1_ERI_O)
      (sa, ia, sa4, ia4, W302aa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W303(j,i,a1,a0) += (    1.00000000) V2(a1,a3,a2,a0) D2(j,a3,i,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) T2(a0,a1,k,a) W303(j,i,a1,a0) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W303aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa1));
  FC_FUNC(g_if_sigma_ooov_ooov_no7_x0_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO7_X0_TYPE1_ERI_O)
    (sa1, ia1, V2_sym.cptr(), W303aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no7_x1_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO7_X1_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), W303aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W309(k,a1,a0,j,i,a2) += (    1.00000000) V2(i,a3,a4,a2) D3(k,a3,a1,j,a0,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) T2(a0,a1,a2,a) W309(k,a1,a0,j,i,a2) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W309aaaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_ooov_ooov_no8_x0_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO8_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W309aaaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no8_x1_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO8_X1_TYPE1_ERI_O)
      (sa, ia, si, ii, T2b.cptr(), W309aaaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W310(j,a1,a0,i) += (    1.00000000) V2(i,a3,a4,a2) D3(j,a1,a4,a2,a3,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) T2(a0,a1,k,a) W310(j,a1,a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W310aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_ooov_ooov_no9_x0_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO9_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W310aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no9_x1_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO9_X1_TYPE1_ERI_O)
      (sa, ia, si, ii, T2b.cptr(), W310aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W311(i,a1,a0,j) += (    1.00000000) V2(j,a3,a4,a2) D3(i,a1,a4,a2,a3,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) T2(a1,a0,k,a) W311(i,a1,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W311aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ooov_ooov_no10_x0_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO10_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W311aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no10_x1_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO10_X1_TYPE1_ERI_O)
      (sa, ia, sj, ij, T2b.cptr(), W311aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W312(k,a1,a0,i,j,a2) += (    1.00000000) V2(j,a3,a4,a2) D3(k,i,a1,a4,a0,a3) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) T2(a1,a0,a2,a) W312(k,a1,a0,i,j,a2) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W312aaaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ooov_ooov_no11_x0_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO11_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W312aaaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no11_x1_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO11_X1_TYPE1_ERI_O)
      (sa, ia, sj, ij, T2b.cptr(), W312aaaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W313(k,a3,a4,a) += (    1.00000000) T2(a1,a2,a0,a) D3(k,a3,a2,a4,a1,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) V2(a4,j,i,a3) W313(k,a3,a4,a) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W313aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa4^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no12_x0_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO12_X0_TYPE1_ERI_O)
      (sa, ia, sa4, ia4, T2b.cptr(), W313aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no12_x1_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO12_X1_TYPE1_ERI_O)
      (sa, ia, sa4, ia4, V2_sym.cptr(), W313aa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W314(a1,a0,j,i) += (    1.00000000) V2(j,a3,i,a2) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) T2(a0,a1,k,a) W314(a1,a0,j,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W314aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ooov_ooov_no13_x0_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO13_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W314aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no13_x1_type1_eri_o,G_IF_SIGMA_OOOV_OOOV_NO13_X1_TYPE1_ERI_O)
      (sa, ia, sj, ij, T2b.cptr(), W314aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@type(2).declaration(begin)
  // --  Title : sigma_ooov_ooov
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W34av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W36av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W42av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W46aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W50aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W60aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W62av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W66aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W68aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W80aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W82av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W84aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W88aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W90aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W92aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W94aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W96av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W100aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W104aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W114aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W116av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W120aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W122aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W134aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W136av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W138aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W142aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W144aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W146aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W148aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W174av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W176av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W182av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W184aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W188aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W192aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W194aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W196aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W198aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W202av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W204aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W210aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W212aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W214aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W216aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W220aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W222av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W224aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W228aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W230aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W232aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W234aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W236av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W238aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W242aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W246aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W248aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W250aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W252aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W256av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W258aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W264aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W266aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W268aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W270aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W274aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W276av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W278aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W282aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W284aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W286aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W288aaav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W11av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W25av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W45av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W49av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W99av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W103av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W151av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W173av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W201av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W207av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W255av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  orz::DTensor W261av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
//-@type(2).declaration(end)

  //*-- Entering to take the type 0 contractions --*//
//-@type(0).contraction(begin)

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W34(a2,v0) += (    1.00000000) T2(a2,a1,v0,a0) D1(a0,a1) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO0_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W34av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W36(a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D1(a0,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO1_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W36av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W42(a3,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a1,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no2_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO2_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W42av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W46(a3,k,a2,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a1,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no3_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO3_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W46aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W50(a3,k,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,k,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no4_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO4_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W50aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W60(j,a4,a2,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a1,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no5_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO5_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W60aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W62(j,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a1,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no6_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO6_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W62av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W66(j,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a2,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no7_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO7_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W66aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W68(a3,i,a0,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a3,a2,i,a1) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no8_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO8_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W68aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W80(i,a4,a2,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a1,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no9_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO9_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W80aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W82(i,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a1,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no10_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO10_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W82av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W84(j,a3,a0,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(j,a1,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no11_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO11_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W84aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W88(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a4,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no12_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO12_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W88aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W90(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a4,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no13_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO13_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W90aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W92(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no14_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO14_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W92aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W94(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no15_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO15_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W94aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W96(a3,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a1,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no16_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO16_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W96av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W100(a3,k,a2,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a1,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no17_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO17_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W100aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W104(a3,k,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,k,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no18_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO18_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W104aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W114(j,a4,a2,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a1,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no19_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO19_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W114aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W116(j,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a1,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no20_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO20_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W116av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W120(j,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a2,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no21_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO21_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W120aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W122(a3,i,a0,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a1,i,a2,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no22_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO22_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W122aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W134(i,a4,a2,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a1,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no23_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO23_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W134aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W136(i,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a1,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no24_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO24_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W136av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W138(j,a3,a0,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a1,j,a2,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no25_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO25_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W138aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W142(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a4,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no26_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO26_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W142aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W144(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a4,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no27_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO27_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W144aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W146(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no28_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO28_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W146aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W148(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no29_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO29_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W148aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W174(a2,v0) += (    1.00000000) T2(a2,a1,v0,a0) D1(a0,a1) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no30_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO30_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W174av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W176(a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D1(a0,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no31_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO31_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W176av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W182(j,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a1,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no32_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO32_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W182av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W184(j,a4,a2,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a1,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no33_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO33_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W184aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W188(j,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a4,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no34_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO34_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W188aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W192(a3,i,a0,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a3,a1,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no35_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO35_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W192aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W194(a3,i,a0,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a3,a2,i,a1) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no36_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO36_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W194aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W196(a3,i,a0,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a3,a1,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no37_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO37_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W196aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W198(a3,i,a0,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a3,a2,i,a1) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no38_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO38_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W198aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W202(a3,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a1,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no39_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO39_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W202av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| W204(a3,k,a2,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a1,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no40_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO40_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W204aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| W210(a3,k,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a2,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no41_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO41_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W210aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| W212(a3,k,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a2,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no42_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO42_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W212aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| W214(a3,k,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,k,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no43_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO43_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W214aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| W216(a3,k,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,k,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no44_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO44_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W216aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| W220(i,a4,a2,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a1,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no45_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO45_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W220aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| W222(i,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a1,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no46_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO46_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W222av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| W224(a3,j,a0,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a3,a1,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no47_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO47_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W224aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   48] -- 
  // |-- [    0] --| W228(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a4,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no48_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO48_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W228aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   49] -- 
  // |-- [    0] --| W230(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a4,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no49_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO49_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W230aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   50] -- 
  // |-- [    0] --| W232(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no50_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO50_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W232aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   51] -- 
  // |-- [    0] --| W234(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) D2(a0,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no51_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO51_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W234aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   52] -- 
  // |-- [    0] --| W236(j,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a1,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no52_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO52_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W236av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   53] -- 
  // |-- [    0] --| W238(j,a4,a2,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a1,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no53_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO53_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W238aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   54] -- 
  // |-- [    0] --| W242(j,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a4,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no54_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO54_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W242aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   55] -- 
  // |-- [    0] --| W246(a3,i,a0,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a1,a3,a2,i) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no55_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO55_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W246aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   56] -- 
  // |-- [    0] --| W248(a3,i,a0,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a1,i,a2,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no56_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO56_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W248aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   57] -- 
  // |-- [    0] --| W250(a3,i,a0,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a1,a3,a2,i) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no57_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO57_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W250aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   58] -- 
  // |-- [    0] --| W252(a3,i,a0,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a1,i,a2,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no58_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO58_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W252aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   59] -- 
  // |-- [    0] --| W256(a3,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a1,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no59_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO59_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W256av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   60] -- 
  // |-- [    0] --| W258(a3,k,a2,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a1,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no60_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO60_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W258aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   61] -- 
  // |-- [    0] --| W264(a3,k,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a2,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no61_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO61_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W264aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   62] -- 
  // |-- [    0] --| W266(a3,k,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a2,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no62_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO62_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W266aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   63] -- 
  // |-- [    0] --| W268(a3,k,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,k,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no63_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO63_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W268aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   64] -- 
  // |-- [    0] --| W270(a3,k,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,k,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no64_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO64_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W270aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   65] -- 
  // |-- [    0] --| W274(i,a4,a2,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a1,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no65_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO65_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W274aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   66] -- 
  // |-- [    0] --| W276(i,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a1,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no66_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO66_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W276av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   67] -- 
  // |-- [    0] --| W278(a3,j,a0,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a1,a3,a2,j) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no67_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO67_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W278aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   68] -- 
  // |-- [    0] --| W282(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a4,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no68_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO68_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W282aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   69] -- 
  // |-- [    0] --| W284(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a4,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no69_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO69_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W284aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   70] -- 
  // |-- [    0] --| W286(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no70_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO70_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W286aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   71] -- 
  // |-- [    0] --| W288(i,a4,a1,v0) += (    1.00000000) T2(a2,a1,v0,a0) C2(a0,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_ooov_no71_x0_type0_eri_v,G_IF_SIGMA_OOOV_OOOV_NO71_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W288aaav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

//-@type(0).contraction(end)

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
  // -- Title : sigma_ooov_ooov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W10(a3,j,a4,v0) += (    1.00000000) T2(a1,a2,a0,v0) D3(a0,a1,a3,a2,j,a4) 
  // |-- [    1] --| W11(j,a) += (    1.00000000) V2(v0,a3,a4,a) W10(a3,j,a4,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W10aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no0_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO0_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W10aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no0_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO0_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W10aaa_sigma_ooov_ooov.cptr(), W11av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W12(a2,a1,a0,a4,a3,a) += (    1.00000000) V2(a,a4,v0,a3) T2(a2,a1,v0,a0) 
  // |-- [    1] --| W13(i,a2,a0,a4,a3,a) += (    1.00000000) D1(i,a1) W12(a2,a1,a0,a4,a3,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) D3(a0,k,a3,a2,j,a4) W13(i,a2,a0,a4,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W12aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    orz::DTensor W13aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO1_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W12aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO1_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W12aaaa_sigma_ooov_ooov.cptr(), W13aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO1_X2_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W13aaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W14(a2,a1,a0,a4,a3,a) += (    1.00000000) V2(a,a4,v0,a3) T2(a2,a1,v0,a0) 
  // |-- [    1] --| W15(i,a1,a0,a4,a3,a) += (    1.00000000) D1(i,a2) W14(a2,a1,a0,a4,a3,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) D3(a0,a1,a3,k,j,a4) W15(i,a1,a0,a4,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W14aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    orz::DTensor W15aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no2_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO2_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W14aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no2_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO2_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W14aaaa_sigma_ooov_ooov.cptr(), W15aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no2_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO2_X2_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W15aaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W17(a3,j,k,v0) += (    1.00000000) T2(a1,a2,a0,v0) D3(a0,a1,a3,a2,j,k) 
  // |-- [    1] --| W16(i,a3,v0,a) += (    1.00000000) V2(v0,a3,a4,a) D1(i,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) W16(i,a3,v0,a) W17(a3,j,k,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W17aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no3_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO3_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W17aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W16aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no3_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO3_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W16aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no3_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO3_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W16aa_sigma_ooov_ooov.cptr(), W17aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W19(a3,i,k,v0) += (    1.00000000) T2(a1,a2,a0,v0) D3(a0,a1,a3,a2,i,k) 
  // |-- [    1] --| W18(j,a3,v0,a) += (    1.00000000) V2(v0,a3,a4,a) D1(j,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (    1.00000000) W18(j,a3,v0,a) W19(a3,i,k,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W19aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no4_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO4_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W19aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W18aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no4_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO4_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W18aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no4_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO4_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W18aa_sigma_ooov_ooov.cptr(), W19aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W20(a2,a1,a0,a4,a3,a) += (    1.00000000) V2(a,a4,v0,a3) T2(a2,a1,v0,a0) 
  // |-- [    1] --| W21(j,a2,a0,a4,a3,a) += (    1.00000000) D1(j,a1) W20(a2,a1,a0,a4,a3,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) D3(a0,a4,a3,a2,i,k) W21(j,a2,a0,a4,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W20aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    orz::DTensor W21aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no5_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO5_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W20aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no5_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO5_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W20aaaa_sigma_ooov_ooov.cptr(), W21aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no5_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO5_X2_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W21aaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W22(a2,a1,a0,a4,a3,a) += (    1.00000000) V2(a,a4,v0,a3) T2(a2,a1,v0,a0) 
  // |-- [    1] --| W23(j,a1,a0,a4,a3,a) += (    1.00000000) D1(j,a2) W22(a2,a1,a0,a4,a3,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) D3(a0,a1,a3,a4,i,k) W23(j,a1,a0,a4,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W22aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    orz::DTensor W23aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no6_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO6_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W22aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no6_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO6_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W22aaaa_sigma_ooov_ooov.cptr(), W23aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no6_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO6_X2_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W23aaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W24(a3,i,a4,v0) += (    1.00000000) T2(a1,a2,a0,v0) D3(a0,a1,a3,a2,i,a4) 
  // |-- [    1] --| W25(i,a) += (    1.00000000) V2(v0,a3,a4,a) W24(a3,i,a4,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W24aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no7_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO7_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W24aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no7_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO7_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W24aaa_sigma_ooov_ooov.cptr(), W25av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W26(a2,a4,v0,a) += (    1.00000000) V2(a,a4,v0,a3) D1(a3,a2) 
  // |-- [    1] --| W27(a1,a0,a4,a) += (    1.00000000) T2(a2,a1,v0,a0) W26(a2,a4,v0,a) 
  // |-- [    2] --| S2(i,j,k,a) += (    1.00000000) D3(a0,a1,j,a4,i,k) W27(a1,a0,a4,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W26aav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no8_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO8_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W26aav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W27aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no8_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO8_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W26aav_sigma_ooov_ooov.cptr(), W27aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no8_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO8_X2_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W27aa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W28(a1,a4,v0,a) += (    1.00000000) V2(a,a4,v0,a3) D1(a3,a1) 
  // |-- [    1] --| W29(a2,a0,a4,a) += (    1.00000000) T2(a2,a1,v0,a0) W28(a1,a4,v0,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) D3(a0,a2,j,a4,i,k) W29(a2,a0,a4,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W28aav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no9_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO9_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W28aav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W29aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no9_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO9_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W28aav_sigma_ooov_ooov.cptr(), W29aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no9_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO9_X2_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W29aa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W31(j,i,k,v0) += (    1.00000000) T2(a1,a2,a0,v0) D3(a0,a1,j,a2,i,k) 
  // |-- [    1] --| W30(v0,a) += (    1.00000000) V2(v0,a3,a4,a) D1(a3,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) W30(v0,a) W31(j,i,k,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W31aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no10_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO10_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W31aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    double W30_sigma_ooov_ooov(0);
    FC_FUNC(g_if_sigma_ooov_ooov_no10_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO10_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), &W30_sigma_ooov_ooov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no10_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO10_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, &W30_sigma_ooov_ooov, W31aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W33(j,i,a4,v0) += (    1.00000000) T2(a1,a2,a0,v0) D3(a0,a1,j,a4,i,a2) 
  // |-- [    1] --| W32(k,a4,v0,a) += (    1.00000000) V2(v0,a3,a4,a) D1(a3,k) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) W32(k,a4,v0,a) W33(j,i,a4,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W33aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no11_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO11_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W33aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W32aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no11_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO11_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W32aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no11_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO11_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W32aa_sigma_ooov_ooov.cptr(), W33aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W35(a4,a3,a2,a) += (    1.00000000) V2(a,a4,v0,a3) W34(a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) D3(a3,a2,j,a4,i,k) W35(a4,a3,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W35aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no12_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO12_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W34av_sigma_ooov_ooov.cptr(), W35aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no12_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO12_X1_TYPE1_ERI_V)
    (sa, ia, W35aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W37(a4,a3,a1,a) += (    1.00000000) V2(a,a4,v0,a3) W36(a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) D3(a3,a1,j,a4,i,k) W37(a4,a3,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W37aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no13_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO13_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W36av_sigma_ooov_ooov.cptr(), W37aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no13_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO13_X1_TYPE1_ERI_V)
    (sa, ia, W37aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W38(a0,a3,v0,a) += (    1.00000000) V2(a,a4,v0,a3) D1(a0,a4) 
  // |-- [    1] --| W39(a2,a1,a3,a) += (    1.00000000) T2(a2,a1,v0,a0) W38(a0,a3,v0,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) D3(a3,a2,j,a1,i,k) W39(a2,a1,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W39aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W38av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no14_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO14_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W38av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no14_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO14_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W38av_sigma_ooov_ooov.cptr(), W39aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  FC_FUNC(g_if_sigma_ooov_ooov_no14_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO14_X2_TYPE1_ERI_V)
    (sa, ia, W39aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W40(a2,a1,a0,a4,a3,a) += (    1.00000000) V2(a,a4,v0,a3) T2(a2,a1,v0,a0) 
  // |-- [    1] --| W41(k,a2,a1,a4,a3,a) += (    1.00000000) D1(a0,k) W40(a2,a1,a0,a4,a3,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) D3(a3,a2,j,a4,i,a1) W41(k,a2,a1,a4,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W41aaaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sa));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W40aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no15_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO15_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W40aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no15_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO15_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W40aaaa_sigma_ooov_ooov.cptr(), W41aaaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  FC_FUNC(g_if_sigma_ooov_ooov_no15_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO15_X2_TYPE1_ERI_V)
    (sa, ia, W41aaaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W43(a4,a) += (    1.00000000) V2(a,a4,v0,a3) W42(a3,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D2(j,a4,i,k) W43(a4,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W43a_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no16_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO16_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W42av_sigma_ooov_ooov.cptr(), W43a_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no16_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO16_X1_TYPE1_ERI_V)
    (sa, ia, W43a_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W44(a3,a4,a2,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,a1,a3,a4) 
  // |-- [    1] --| W45(a2,a) += (    1.00000000) V2(v0,a3,a4,a) W44(a3,a4,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W44aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no17_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO17_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W44aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no17_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO17_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W44aaa_sigma_ooov_ooov.cptr(), W45av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W47(a4,k,a2,a) += (    1.00000000) V2(a,a4,v0,a3) W46(a3,k,a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(j,a4,i,a2) W47(a4,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W47aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no18_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO18_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W46aaav_sigma_ooov_ooov.cptr(), W47aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no18_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO18_X1_TYPE1_ERI_V)
    (sa, ia, W47aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W48(a3,a4,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,a4,a3,a2) 
  // |-- [    1] --| W49(a1,a) += (    1.00000000) V2(v0,a3,a4,a) W48(a3,a4,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W48aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no19_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO19_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W48aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no19_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO19_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W48aaa_sigma_ooov_ooov.cptr(), W49av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W51(a4,k,a1,a) += (    1.00000000) V2(a,a4,v0,a3) W50(a3,k,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(j,a4,i,a1) W51(a4,k,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W51aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no20_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO20_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W50aaav_sigma_ooov_ooov.cptr(), W51aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no20_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO20_X1_TYPE1_ERI_V)
    (sa, ia, W51aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W53(j,i,a0,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(j,a1,i,a2) 
  // |-- [    1] --| W52(a0,k,v0,a) += (    1.00000000) V2(v0,a3,a4,a) D2(a0,a4,a3,k) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.33333333) W52(a0,k,v0,a) W53(j,i,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W53aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no21_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO21_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W53aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W52aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no21_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO21_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W52aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no21_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO21_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W52aa_sigma_ooov_ooov.cptr(), W53aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W55(j,i,a0,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(j,a2,i,a1) 
  // |-- [    1] --| W54(a0,k,v0,a) += (    1.00000000) V2(v0,a3,a4,a) D2(a0,a4,a3,k) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.16666667) W54(a0,k,v0,a) W55(j,i,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W55aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no22_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO22_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W55aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W54aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no22_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO22_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W54aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no22_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO22_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W54aa_sigma_ooov_ooov.cptr(), W55aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W57(j,i,a0,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(j,a1,i,a2) 
  // |-- [    1] --| W56(a0,k,v0,a) += (    1.00000000) V2(v0,a3,a4,a) D2(a0,k,a3,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.16666667) W56(a0,k,v0,a) W57(j,i,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W57aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no23_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO23_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W57aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W56aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no23_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO23_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W56aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no23_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO23_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W56aa_sigma_ooov_ooov.cptr(), W57aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W59(j,i,a0,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(j,a2,i,a1) 
  // |-- [    1] --| W58(a0,k,v0,a) += (    1.00000000) V2(v0,a3,a4,a) D2(a0,k,a3,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.33333333) W58(a0,k,v0,a) W59(j,i,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W59aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no24_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO24_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W59aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W58aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no24_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO24_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W58aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no24_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO24_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W58aa_sigma_ooov_ooov.cptr(), W59aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W61(a3,j,a2,a) += (    1.00000000) V2(a,a4,v0,a3) W60(j,a4,a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D2(a3,a2,i,k) W61(a3,j,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W61aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no25_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO25_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W60aaav_sigma_ooov_ooov.cptr(), W61aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no25_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO25_X1_TYPE1_ERI_V)
    (sa, ia, W61aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W63(a4,a3,j,a) += (    1.00000000) V2(a,a4,v0,a3) W62(j,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(a3,a4,i,k) W63(a4,a3,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W63aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no26_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO26_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W62av_sigma_ooov_ooov.cptr(), W63aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no26_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO26_X1_TYPE1_ERI_V)
    (sa, ia, W63aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W65(j,k,a2,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,a1,j,k) 
  // |-- [    1] --| W64(i,a2,v0,a) += (    1.00000000) V2(v0,a3,a4,a) D2(a3,a2,i,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.50000000) W64(i,a2,v0,a) W65(j,k,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W65aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no27_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO27_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W65aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W64aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no27_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO27_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W64aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no27_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO27_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W64aa_sigma_ooov_ooov.cptr(), W65aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W67(a3,j,a1,a) += (    1.00000000) V2(a,a4,v0,a3) W66(j,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(a3,a1,i,k) W67(a3,j,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W67aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no28_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO28_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W66aaav_sigma_ooov_ooov.cptr(), W67aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no28_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO28_X1_TYPE1_ERI_V)
    (sa, ia, W67aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W69(a4,i,a0,a) += (    1.00000000) V2(a,a4,v0,a3) W68(a3,i,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(a0,k,j,a4) W69(a4,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W69aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no29_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO29_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W68aaav_sigma_ooov_ooov.cptr(), W69aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no29_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO29_X1_TYPE1_ERI_V)
    (sa, ia, W69aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W71(j,k,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,a2,j,k) 
  // |-- [    1] --| W70(i,a1,v0,a) += (    1.00000000) V2(v0,a3,a4,a) D2(a3,a1,i,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.33333333) W70(i,a1,v0,a) W71(j,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W71aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no30_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO30_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W71aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W70aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no30_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO30_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W70aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no30_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO30_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W70aa_sigma_ooov_ooov.cptr(), W71aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W73(j,k,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,a2,j,k) 
  // |-- [    1] --| W72(i,a1,v0,a) += (    1.00000000) V2(v0,a3,a4,a) D2(a3,a4,i,a1) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.16666667) W72(i,a1,v0,a) W73(j,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W73aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no31_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO31_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W73aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W72aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no31_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO31_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W72aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no31_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO31_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W72aa_sigma_ooov_ooov.cptr(), W73aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W75(j,k,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,k,j,a2) 
  // |-- [    1] --| W74(i,a1,v0,a) += (    1.00000000) V2(v0,a3,a4,a) D2(a3,a1,i,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.16666667) W74(i,a1,v0,a) W75(j,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W75aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no32_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO32_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W75aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W74aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no32_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO32_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W74aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no32_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO32_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W74aa_sigma_ooov_ooov.cptr(), W75aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W77(j,k,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,k,j,a2) 
  // |-- [    1] --| W76(i,a1,v0,a) += (    1.00000000) V2(v0,a3,a4,a) D2(a3,a4,i,a1) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.33333333) W76(i,a1,v0,a) W77(j,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W77aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no33_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO33_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W77aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W76aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no33_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO33_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W76aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no33_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO33_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W76aa_sigma_ooov_ooov.cptr(), W77aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W79(i,k,a2,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,a1,i,k) 
  // |-- [    1] --| W78(j,a2,v0,a) += (    1.00000000) V2(v0,a3,a4,a) D2(j,a4,a3,a2) 
  // |-- [    2] --| S2(i,j,k,a) += (   -1.00000000) W78(j,a2,v0,a) W79(i,k,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W79aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no34_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO34_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W79aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W78aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no34_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO34_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W78aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no34_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO34_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W78aa_sigma_ooov_ooov.cptr(), W79aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W81(a3,i,a2,a) += (    1.00000000) V2(a,a4,v0,a3) W80(i,a4,a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(j,k,a3,a2) W81(a3,i,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W81aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no35_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO35_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W80aaav_sigma_ooov_ooov.cptr(), W81aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no35_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO35_X1_TYPE1_ERI_V)
    (sa, ia, W81aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W83(a4,a3,i,a) += (    1.00000000) V2(a,a4,v0,a3) W82(i,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(j,a4,a3,k) W83(a4,a3,i,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W83aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no36_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO36_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W82av_sigma_ooov_ooov.cptr(), W83aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no36_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO36_X1_TYPE1_ERI_V)
    (sa, ia, W83aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W85(a4,j,a0,a) += (    1.00000000) V2(a,a4,v0,a3) W84(j,a3,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(a0,a4,i,k) W85(a4,j,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W85aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no37_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO37_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W84aaav_sigma_ooov_ooov.cptr(), W85aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no37_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO37_X1_TYPE1_ERI_V)
    (sa, ia, W85aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W87(i,k,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,a2,i,k) 
  // |-- [    1] --| W86(j,a1,v0,a) += (    1.00000000) V2(v0,a3,a4,a) D2(j,a4,a3,a1) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.50000000) W86(j,a1,v0,a) W87(i,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W87aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no38_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO38_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W87aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W86aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no38_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO38_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W86aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no38_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO38_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W86aa_sigma_ooov_ooov.cptr(), W87aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W89(a3,i,a1,a) += (    1.00000000) V2(a,a4,v0,a3) W88(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.33333333) D2(j,a1,a3,k) W89(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W89aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no39_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO39_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W88aaav_sigma_ooov_ooov.cptr(), W89aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no39_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO39_X1_TYPE1_ERI_V)
    (sa, ia, W89aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| W91(a3,i,a1,a) += (    1.00000000) V2(a,a4,v0,a3) W90(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.16666667) D2(j,k,a3,a1) W91(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W91aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no40_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO40_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W90aaav_sigma_ooov_ooov.cptr(), W91aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no40_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO40_X1_TYPE1_ERI_V)
    (sa, ia, W91aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| W93(a3,i,a1,a) += (    1.00000000) V2(a,a4,v0,a3) W92(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.16666667) D2(j,a1,a3,k) W93(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W93aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no41_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO41_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W92aaav_sigma_ooov_ooov.cptr(), W93aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no41_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO41_X1_TYPE1_ERI_V)
    (sa, ia, W93aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| W95(a3,i,a1,a) += (    1.00000000) V2(a,a4,v0,a3) W94(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.33333333) D2(j,k,a3,a1) W95(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W95aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no42_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO42_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W94aaav_sigma_ooov_ooov.cptr(), W95aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no42_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO42_X1_TYPE1_ERI_V)
    (sa, ia, W95aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| W97(a4,a) += (    1.00000000) V2(a,a4,v0,a3) W96(a3,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    2.00000000) C2(a4,j,k,i) W97(a4,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W97a_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no43_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO43_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W96av_sigma_ooov_ooov.cptr(), W97a_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no43_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO43_X1_TYPE1_ERI_V)
    (sa, ia, W97a_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| W98(a3,a4,a2,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,a1,a3,a4) 
  // |-- [    1] --| W99(a2,a) += (    1.00000000) V2(v0,a3,a4,a) W98(a3,a4,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W98aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no44_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO44_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W98aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no44_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO44_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W98aaa_sigma_ooov_ooov.cptr(), W99av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| W101(a4,k,a2,a) += (    1.00000000) V2(a,a4,v0,a3) W100(a3,k,a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a2,i,a4,j) W101(a4,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W101aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no45_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO45_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W100aaav_sigma_ooov_ooov.cptr(), W101aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no45_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO45_X1_TYPE1_ERI_V)
    (sa, ia, W101aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| W102(a3,a4,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,a4,a3,a2) 
  // |-- [    1] --| W103(a1,a) += (    1.00000000) V2(v0,a3,a4,a) W102(a3,a4,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W102aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no46_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO46_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W102aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no46_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO46_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W102aaa_sigma_ooov_ooov.cptr(), W103av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| W105(a4,k,a1,a) += (    1.00000000) V2(a,a4,v0,a3) W104(a3,k,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a1,i,a4,j) W105(a4,k,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W105aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no47_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO47_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W104aaav_sigma_ooov_ooov.cptr(), W105aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no47_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO47_X1_TYPE1_ERI_V)
    (sa, ia, W105aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   48] -- 
  // |-- [    0] --| W107(j,i,a0,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a1,j,a2,i) 
  // |-- [    1] --| W106(a0,k,v0,a) += (    1.00000000) V2(v0,a3,a4,a) C2(a0,a4,a3,k) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.66666667) W106(a0,k,v0,a) W107(j,i,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W107aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no48_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO48_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W107aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W106aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no48_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO48_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W106aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no48_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO48_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W106aa_sigma_ooov_ooov.cptr(), W107aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   49] -- 
  // |-- [    0] --| W109(j,i,a0,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a1,i,a2,j) 
  // |-- [    1] --| W108(a0,k,v0,a) += (    1.00000000) V2(v0,a3,a4,a) C2(a0,a4,a3,k) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.33333333) W108(a0,k,v0,a) W109(j,i,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W109aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no49_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO49_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W109aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W108aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no49_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO49_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W108aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no49_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO49_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W108aa_sigma_ooov_ooov.cptr(), W109aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   50] -- 
  // |-- [    0] --| W111(j,i,a0,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a1,j,a2,i) 
  // |-- [    1] --| W110(a0,k,v0,a) += (    1.00000000) V2(v0,a3,a4,a) C2(a0,k,a3,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.33333333) W110(a0,k,v0,a) W111(j,i,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W111aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no50_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO50_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W111aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W110aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no50_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO50_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W110aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no50_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO50_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W110aa_sigma_ooov_ooov.cptr(), W111aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   51] -- 
  // |-- [    0] --| W113(j,i,a0,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a1,i,a2,j) 
  // |-- [    1] --| W112(a0,k,v0,a) += (    1.00000000) V2(v0,a3,a4,a) C2(a0,k,a3,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.66666667) W112(a0,k,v0,a) W113(j,i,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W113aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no51_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO51_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W113aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W112aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no51_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO51_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W112aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no51_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO51_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W112aa_sigma_ooov_ooov.cptr(), W113aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   52] -- 
  // |-- [    0] --| W115(a3,j,a2,a) += (    1.00000000) V2(a,a4,v0,a3) W114(j,a4,a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    2.00000000) C2(a2,a3,k,i) W115(a3,j,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W115aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no52_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO52_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W114aaav_sigma_ooov_ooov.cptr(), W115aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no52_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO52_X1_TYPE1_ERI_V)
    (sa, ia, W115aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   53] -- 
  // |-- [    0] --| W117(a4,a3,j,a) += (    1.00000000) V2(a,a4,v0,a3) W116(j,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a3,a4,i,k) W117(a4,a3,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W117aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no53_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO53_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W116av_sigma_ooov_ooov.cptr(), W117aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no53_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO53_X1_TYPE1_ERI_V)
    (sa, ia, W117aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   54] -- 
  // |-- [    0] --| W119(j,k,a2,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,a1,j,k) 
  // |-- [    1] --| W118(i,a2,v0,a) += (    1.00000000) V2(v0,a3,a4,a) C2(a2,a3,a4,i) 
  // |-- [    2] --| S2(i,j,k,a) += (   -1.00000000) W118(i,a2,v0,a) W119(j,k,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W119aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no54_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO54_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W119aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W118aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no54_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO54_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W118aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no54_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO54_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W118aa_sigma_ooov_ooov.cptr(), W119aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   55] -- 
  // |-- [    0] --| W121(a3,j,a1,a) += (    1.00000000) V2(a,a4,v0,a3) W120(j,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a1,a3,k,i) W121(a3,j,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W121aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no55_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO55_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W120aaav_sigma_ooov_ooov.cptr(), W121aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no55_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO55_X1_TYPE1_ERI_V)
    (sa, ia, W121aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   56] -- 
  // |-- [    0] --| W123(a4,i,a0,a) += (    1.00000000) V2(a,a4,v0,a3) W122(a3,i,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a0,k,j,a4) W123(a4,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W123aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no56_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO56_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W122aaav_sigma_ooov_ooov.cptr(), W123aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no56_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO56_X1_TYPE1_ERI_V)
    (sa, ia, W123aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   57] -- 
  // |-- [    0] --| W125(j,k,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,a2,j,k) 
  // |-- [    1] --| W124(i,a1,v0,a) += (    1.00000000) V2(v0,a3,a4,a) C2(a1,a3,a4,i) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.66666667) W124(i,a1,v0,a) W125(j,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W125aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no57_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO57_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W125aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W124aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no57_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO57_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W124aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no57_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO57_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W124aa_sigma_ooov_ooov.cptr(), W125aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   58] -- 
  // |-- [    0] --| W127(j,k,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,a2,j,k) 
  // |-- [    1] --| W126(i,a1,v0,a) += (    1.00000000) V2(v0,a3,a4,a) C2(a1,i,a4,a3) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.33333333) W126(i,a1,v0,a) W127(j,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W127aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no58_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO58_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W127aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W126aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no58_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO58_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W126aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no58_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO58_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W126aa_sigma_ooov_ooov.cptr(), W127aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   59] -- 
  // |-- [    0] --| W129(j,k,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,k,j,a2) 
  // |-- [    1] --| W128(i,a1,v0,a) += (    1.00000000) V2(v0,a3,a4,a) C2(a1,a3,a4,i) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.33333333) W128(i,a1,v0,a) W129(j,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W129aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no59_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO59_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W129aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W128aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no59_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO59_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W128aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no59_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO59_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W128aa_sigma_ooov_ooov.cptr(), W129aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   60] -- 
  // |-- [    0] --| W131(j,k,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,k,j,a2) 
  // |-- [    1] --| W130(i,a1,v0,a) += (    1.00000000) V2(v0,a3,a4,a) C2(a1,i,a4,a3) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.66666667) W130(i,a1,v0,a) W131(j,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W131aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no60_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO60_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W131aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W130aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no60_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO60_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W130aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no60_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO60_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W130aa_sigma_ooov_ooov.cptr(), W131aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   61] -- 
  // |-- [    0] --| W133(i,k,a2,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,a1,i,k) 
  // |-- [    1] --| W132(j,a2,v0,a) += (    1.00000000) V2(v0,a3,a4,a) C2(a2,a3,a4,j) 
  // |-- [    2] --| S2(i,j,k,a) += (    2.00000000) W132(j,a2,v0,a) W133(i,k,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W133aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no61_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO61_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W133aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W132aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no61_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO61_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W132aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no61_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO61_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W132aa_sigma_ooov_ooov.cptr(), W133aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   62] -- 
  // |-- [    0] --| W135(a3,i,a2,a) += (    1.00000000) V2(a,a4,v0,a3) W134(i,a4,a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a2,a3,k,j) W135(a3,i,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W135aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no62_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO62_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W134aaav_sigma_ooov_ooov.cptr(), W135aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no62_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO62_X1_TYPE1_ERI_V)
    (sa, ia, W135aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   63] -- 
  // |-- [    0] --| W137(a4,a3,i,a) += (    1.00000000) V2(a,a4,v0,a3) W136(i,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a3,k,j,a4) W137(a4,a3,i,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W137aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no63_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO63_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W136av_sigma_ooov_ooov.cptr(), W137aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no63_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO63_X1_TYPE1_ERI_V)
    (sa, ia, W137aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   64] -- 
  // |-- [    0] --| W139(a4,j,a0,a) += (    1.00000000) V2(a,a4,v0,a3) W138(j,a3,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a0,a4,i,k) W139(a4,j,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W139aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no64_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO64_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W138aaav_sigma_ooov_ooov.cptr(), W139aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no64_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO64_X1_TYPE1_ERI_V)
    (sa, ia, W139aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   65] -- 
  // |-- [    0] --| W141(i,k,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,a2,i,k) 
  // |-- [    1] --| W140(j,a1,v0,a) += (    1.00000000) V2(v0,a3,a4,a) C2(a1,a3,a4,j) 
  // |-- [    2] --| S2(i,j,k,a) += (   -1.00000000) W140(j,a1,v0,a) W141(i,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W141aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no65_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO65_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W141aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W140aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no65_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO65_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W140aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no65_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO65_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W140aa_sigma_ooov_ooov.cptr(), W141aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   66] -- 
  // |-- [    0] --| W143(a3,i,a1,a) += (    1.00000000) V2(a,a4,v0,a3) W142(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.66666667) C2(a1,j,k,a3) W143(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W143aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no66_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO66_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W142aaav_sigma_ooov_ooov.cptr(), W143aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no66_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO66_X1_TYPE1_ERI_V)
    (sa, ia, W143aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   67] -- 
  // |-- [    0] --| W145(a3,i,a1,a) += (    1.00000000) V2(a,a4,v0,a3) W144(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.33333333) C2(a1,a3,k,j) W145(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W145aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no67_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO67_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W144aaav_sigma_ooov_ooov.cptr(), W145aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no67_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO67_X1_TYPE1_ERI_V)
    (sa, ia, W145aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   68] -- 
  // |-- [    0] --| W147(a3,i,a1,a) += (    1.00000000) V2(a,a4,v0,a3) W146(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.33333333) C2(a1,j,k,a3) W147(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W147aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no68_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO68_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W146aaav_sigma_ooov_ooov.cptr(), W147aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no68_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO68_X1_TYPE1_ERI_V)
    (sa, ia, W147aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   69] -- 
  // |-- [    0] --| W149(a3,i,a1,a) += (    1.00000000) V2(a,a4,v0,a3) W148(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.66666667) C2(a1,a3,k,j) W149(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W149aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no69_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO69_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W148aaav_sigma_ooov_ooov.cptr(), W149aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no69_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO69_X1_TYPE1_ERI_V)
    (sa, ia, W149aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   70] -- 
  // |-- [    0] --| W150(j,a3,a4,v0) += (    1.00000000) T2(a1,a2,a0,v0) D3(a0,a1,j,a2,a3,a4) 
  // |-- [    1] --| W151(j,a) += (    1.00000000) V2(v0,a,a4,a3) W150(j,a3,a4,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W150aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no70_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO70_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W150aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no70_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO70_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W150aaa_sigma_ooov_ooov.cptr(), W151av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   71] -- 
  // |-- [    0] --| W152(a2,a1,a0,a4,a3,a) += (    1.00000000) V2(a,v0,a4,a3) T2(a2,a1,v0,a0) 
  // |-- [    1] --| W153(i,a2,a0,a4,a3,a) += (    1.00000000) D1(i,a1) W152(a2,a1,a0,a4,a3,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) D3(a0,k,j,a2,a3,a4) W153(i,a2,a0,a4,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W152aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    orz::DTensor W153aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no71_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO71_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W152aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no71_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO71_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W152aaaa_sigma_ooov_ooov.cptr(), W153aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no71_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO71_X2_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W153aaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   72] -- 
  // |-- [    0] --| W154(a2,a1,a0,a4,a3,a) += (    1.00000000) V2(a,v0,a4,a3) T2(a2,a1,v0,a0) 
  // |-- [    1] --| W155(i,a1,a0,a4,a3,a) += (    1.00000000) D1(i,a2) W154(a2,a1,a0,a4,a3,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) D3(a0,a1,j,k,a3,a4) W155(i,a1,a0,a4,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W154aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    orz::DTensor W155aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no72_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO72_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W154aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no72_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO72_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W154aaaa_sigma_ooov_ooov.cptr(), W155aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no72_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO72_X2_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W155aaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   73] -- 
  // |-- [    0] --| W157(j,a3,k,v0) += (    1.00000000) T2(a1,a2,a0,v0) D3(a0,a1,j,a2,a3,k) 
  // |-- [    1] --| W156(i,a3,v0,a) += (    1.00000000) V2(v0,a,a4,a3) D1(i,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) W156(i,a3,v0,a) W157(j,a3,k,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W157aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no73_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO73_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W157aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W156aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no73_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO73_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W156aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no73_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO73_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W156aa_sigma_ooov_ooov.cptr(), W157aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   74] -- 
  // |-- [    0] --| W159(j,i,k,v0) += (    1.00000000) T2(a1,a2,a0,v0) D3(a0,a1,j,a2,i,k) 
  // |-- [    1] --| W158(v0,a) += (    1.00000000) V2(v0,a,a4,a3) D1(a3,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (    1.00000000) W158(v0,a) W159(j,i,k,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W159aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no74_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO74_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W159aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    double W158_sigma_ooov_ooov(0);
    FC_FUNC(g_if_sigma_ooov_ooov_no74_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO74_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), &W158_sigma_ooov_ooov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no74_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO74_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, &W158_sigma_ooov_ooov, W159aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   75] -- 
  // |-- [    0] --| W160(a1,a4,v0,a) += (    1.00000000) V2(a,v0,a4,a3) D1(a3,a1) 
  // |-- [    1] --| W161(a2,a0,a4,a) += (    1.00000000) T2(a2,a1,v0,a0) W160(a1,a4,v0,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) D3(a0,a4,j,a2,i,k) W161(a2,a0,a4,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W160aav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no75_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO75_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W160aav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W161aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no75_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO75_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W160aav_sigma_ooov_ooov.cptr(), W161aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no75_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO75_X2_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W161aa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   76] -- 
  // |-- [    0] --| W162(a2,a4,v0,a) += (    1.00000000) V2(a,v0,a4,a3) D1(a3,a2) 
  // |-- [    1] --| W163(a1,a0,a4,a) += (    1.00000000) T2(a2,a1,v0,a0) W162(a2,a4,v0,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) D3(a0,a1,j,a4,i,k) W163(a1,a0,a4,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W162aav_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no76_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO76_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W162aav_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W163aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no76_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO76_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W162aav_sigma_ooov_ooov.cptr(), W163aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no76_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO76_X2_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W163aa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   77] -- 
  // |-- [    0] --| W165(j,i,a4,v0) += (    1.00000000) T2(a1,a2,a0,v0) D3(a0,a1,j,a2,i,a4) 
  // |-- [    1] --| W164(k,a4,v0,a) += (    1.00000000) V2(v0,a,a4,a3) D1(a3,k) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) W164(k,a4,v0,a) W165(j,i,a4,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W165aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no77_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO77_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W165aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W164aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no77_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO77_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W164aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no77_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO77_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W164aa_sigma_ooov_ooov.cptr(), W165aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   78] -- 
  // |-- [    0] --| W166(a2,a1,a0,a4,a3,a) += (    1.00000000) V2(a,v0,a4,a3) T2(a2,a1,v0,a0) 
  // |-- [    1] --| W167(j,a1,a0,a4,a3,a) += (    1.00000000) D1(j,a2) W166(a2,a1,a0,a4,a3,a) 
  // |-- [    2] --| S2(i,j,k,a) += (    1.00000000) D3(a0,a1,a3,a4,i,k) W167(j,a1,a0,a4,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W166aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    orz::DTensor W167aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no78_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO78_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W166aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no78_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO78_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W166aaaa_sigma_ooov_ooov.cptr(), W167aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no78_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO78_X2_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W167aaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   79] -- 
  // |-- [    0] --| W168(a2,a1,a0,a4,a3,a) += (    1.00000000) V2(a,v0,a4,a3) T2(a2,a1,v0,a0) 
  // |-- [    1] --| W169(j,a2,a0,a4,a3,a) += (    1.00000000) D1(j,a1) W168(a2,a1,a0,a4,a3,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) D3(a0,a2,a3,a4,i,k) W169(j,a2,a0,a4,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W168aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    orz::DTensor W169aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no79_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO79_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W168aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no79_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO79_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W168aaaa_sigma_ooov_ooov.cptr(), W169aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no79_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO79_X2_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W169aaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   80] -- 
  // |-- [    0] --| W171(a3,i,k,v0) += (    1.00000000) T2(a1,a2,a0,v0) D3(a0,a1,a3,a2,i,k) 
  // |-- [    1] --| W170(j,a3,v0,a) += (    1.00000000) V2(v0,a,a4,a3) D1(j,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) W170(j,a3,v0,a) W171(a3,i,k,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W171aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no80_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO80_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W171aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W170aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no80_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO80_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W170aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no80_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO80_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W170aa_sigma_ooov_ooov.cptr(), W171aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   81] -- 
  // |-- [    0] --| W172(a3,i,a4,v0) += (    1.00000000) T2(a1,a2,a0,v0) D3(a0,a1,a3,a4,i,a2) 
  // |-- [    1] --| W173(i,a) += (    1.00000000) V2(v0,a,a4,a3) W172(a3,i,a4,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W172aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no81_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO81_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W172aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no81_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO81_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W172aaa_sigma_ooov_ooov.cptr(), W173av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   82] -- 
  // |-- [    0] --| W175(a4,a3,a2,a) += (    1.00000000) V2(a,v0,a4,a3) W174(a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) D3(j,a2,a3,a4,i,k) W175(a4,a3,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W175aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no82_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO82_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W174av_sigma_ooov_ooov.cptr(), W175aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no82_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO82_X1_TYPE1_ERI_V)
    (sa, ia, W175aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   83] -- 
  // |-- [    0] --| W177(a4,a3,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W176(a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) D3(j,a1,a3,a4,i,k) W177(a4,a3,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W177aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no83_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO83_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W176av_sigma_ooov_ooov.cptr(), W177aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no83_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO83_X1_TYPE1_ERI_V)
    (sa, ia, W177aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   84] -- 
  // |-- [    0] --| W178(a0,a3,v0,a) += (    1.00000000) V2(a,v0,a4,a3) D1(a0,a4) 
  // |-- [    1] --| W179(a2,a1,a3,a) += (    1.00000000) T2(a2,a1,v0,a0) W178(a0,a3,v0,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) D3(j,a2,a3,a1,i,k) W179(a2,a1,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W179aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W178av_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no84_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO84_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W178av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no84_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO84_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W178av_sigma_ooov_ooov.cptr(), W179aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  FC_FUNC(g_if_sigma_ooov_ooov_no84_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO84_X2_TYPE1_ERI_V)
    (sa, ia, W179aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   85] -- 
  // |-- [    0] --| W180(a2,a1,a0,a4,a3,a) += (    1.00000000) V2(a,v0,a4,a3) T2(a2,a1,v0,a0) 
  // |-- [    1] --| W181(k,a2,a1,a4,a3,a) += (    1.00000000) D1(a0,k) W180(a2,a1,a0,a4,a3,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -0.50000000) D3(j,a2,a3,a4,i,a1) W181(k,a2,a1,a4,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W181aaaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sa));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W180aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no85_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO85_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W180aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no85_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO85_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W180aaaa_sigma_ooov_ooov.cptr(), W181aaaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  FC_FUNC(g_if_sigma_ooov_ooov_no85_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO85_X2_TYPE1_ERI_V)
    (sa, ia, W181aaaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   86] -- 
  // |-- [    0] --| W183(a4,a3,j,a) += (    1.00000000) V2(a,v0,a4,a3) W182(j,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D2(a3,a4,i,k) W183(a4,a3,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W183aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no86_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO86_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W182av_sigma_ooov_ooov.cptr(), W183aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no86_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO86_X1_TYPE1_ERI_V)
    (sa, ia, W183aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   87] -- 
  // |-- [    0] --| W185(a3,j,a2,a) += (    1.00000000) V2(a,v0,a4,a3) W184(j,a4,a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(a3,a2,i,k) W185(a3,j,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W185aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no87_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO87_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W184aaav_sigma_ooov_ooov.cptr(), W185aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no87_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO87_X1_TYPE1_ERI_V)
    (sa, ia, W185aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   88] -- 
  // |-- [    0] --| W187(j,k,a2,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,a1,j,k) 
  // |-- [    1] --| W186(i,a2,v0,a) += (    1.00000000) V2(v0,a,a4,a3) D2(a3,a4,i,a2) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.50000000) W186(i,a2,v0,a) W187(j,k,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W187aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no88_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO88_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W187aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W186aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no88_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO88_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W186aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no88_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO88_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W186aa_sigma_ooov_ooov.cptr(), W187aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   89] -- 
  // |-- [    0] --| W189(a3,j,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W188(j,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(a3,a1,i,k) W189(a3,j,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W189aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no89_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO89_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W188aaav_sigma_ooov_ooov.cptr(), W189aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no89_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO89_X1_TYPE1_ERI_V)
    (sa, ia, W189aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   90] -- 
  // |-- [    0] --| W191(j,k,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,k,j,a2) 
  // |-- [    1] --| W190(i,a1,v0,a) += (    1.00000000) V2(v0,a,a4,a3) D2(a3,a4,i,a1) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.50000000) W190(i,a1,v0,a) W191(j,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W191aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no90_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO90_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W191aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W190aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no90_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO90_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W190aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no90_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO90_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W190aa_sigma_ooov_ooov.cptr(), W191aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   91] -- 
  // |-- [    0] --| W193(a4,i,a0,a) += (    1.00000000) V2(a,v0,a4,a3) W192(a3,i,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.33333333) D2(a0,a4,j,k) W193(a4,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W193aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no91_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO91_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W192aaav_sigma_ooov_ooov.cptr(), W193aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no91_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO91_X1_TYPE1_ERI_V)
    (sa, ia, W193aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   92] -- 
  // |-- [    0] --| W195(a4,i,a0,a) += (    1.00000000) V2(a,v0,a4,a3) W194(a3,i,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.16666667) D2(a0,a4,j,k) W195(a4,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W195aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no92_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO92_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W194aaav_sigma_ooov_ooov.cptr(), W195aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no92_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO92_X1_TYPE1_ERI_V)
    (sa, ia, W195aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   93] -- 
  // |-- [    0] --| W197(a4,i,a0,a) += (    1.00000000) V2(a,v0,a4,a3) W196(a3,i,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.16666667) D2(a0,k,j,a4) W197(a4,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W197aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no93_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO93_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W196aaav_sigma_ooov_ooov.cptr(), W197aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no93_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO93_X1_TYPE1_ERI_V)
    (sa, ia, W197aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   94] -- 
  // |-- [    0] --| W199(a4,i,a0,a) += (    1.00000000) V2(a,v0,a4,a3) W198(a3,i,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.33333333) D2(a0,k,j,a4) W199(a4,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W199aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no94_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO94_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W198aaav_sigma_ooov_ooov.cptr(), W199aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no94_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO94_X1_TYPE1_ERI_V)
    (sa, ia, W199aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   95] -- 
  // |-- [    0] --| W200(a3,a4,a2,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,a1,a3,a4) 
  // |-- [    1] --| W201(a2,a) += (    1.00000000) V2(v0,a,a4,a3) W200(a3,a4,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W200aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no95_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO95_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W200aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no95_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO95_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W200aaa_sigma_ooov_ooov.cptr(), W201av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   96] -- 
  // |-- [    0] --| W203(a4,a) += (    1.00000000) V2(a,v0,a4,a3) W202(a3,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(j,a4,i,k) W203(a4,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W203a_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no96_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO96_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W202av_sigma_ooov_ooov.cptr(), W203a_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no96_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO96_X1_TYPE1_ERI_V)
    (sa, ia, W203a_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   97] -- 
  // |-- [    0] --| W205(a4,k,a2,a) += (    1.00000000) V2(a,v0,a4,a3) W204(a3,k,a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(j,a2,i,a4) W205(a4,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W205aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no97_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO97_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W204aaav_sigma_ooov_ooov.cptr(), W205aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no97_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO97_X1_TYPE1_ERI_V)
    (sa, ia, W205aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   98] -- 
  // |-- [    0] --| W206(a3,a4,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,a2,a3,a4) 
  // |-- [    1] --| W207(a1,a) += (    1.00000000) V2(v0,a,a4,a3) W206(a3,a4,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W206aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no98_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO98_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W206aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no98_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO98_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W206aaa_sigma_ooov_ooov.cptr(), W207av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   99] -- 
  // |-- [    0] --| W209(j,i,a0,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(j,a2,i,a1) 
  // |-- [    1] --| W208(a0,k,v0,a) += (    1.00000000) V2(v0,a,a4,a3) D2(a0,k,a3,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.50000000) W208(a0,k,v0,a) W209(j,i,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W209aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no99_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO99_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W209aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W208aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no99_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO99_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W208aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no99_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO99_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W208aa_sigma_ooov_ooov.cptr(), W209aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  100] -- 
  // |-- [    0] --| W211(a4,k,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W210(a3,k,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.33333333) D2(j,a1,i,a4) W211(a4,k,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W211aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no100_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO100_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W210aaav_sigma_ooov_ooov.cptr(), W211aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no100_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO100_X1_TYPE1_ERI_V)
    (sa, ia, W211aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  101] -- 
  // |-- [    0] --| W213(a4,k,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W212(a3,k,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.16666667) D2(j,a4,i,a1) W213(a4,k,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W213aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no101_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO101_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W212aaav_sigma_ooov_ooov.cptr(), W213aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no101_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO101_X1_TYPE1_ERI_V)
    (sa, ia, W213aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  102] -- 
  // |-- [    0] --| W215(a4,k,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W214(a3,k,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.16666667) D2(j,a1,i,a4) W215(a4,k,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W215aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no102_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO102_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W214aaav_sigma_ooov_ooov.cptr(), W215aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no102_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO102_X1_TYPE1_ERI_V)
    (sa, ia, W215aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  103] -- 
  // |-- [    0] --| W217(a4,k,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W216(a3,k,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.33333333) D2(j,a4,i,a1) W217(a4,k,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W217aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no103_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO103_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W216aaav_sigma_ooov_ooov.cptr(), W217aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no103_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO103_X1_TYPE1_ERI_V)
    (sa, ia, W217aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  104] -- 
  // |-- [    0] --| W219(i,k,a2,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,a1,i,k) 
  // |-- [    1] --| W218(j,a2,v0,a) += (    1.00000000) V2(v0,a,a4,a3) D2(a3,a4,j,a2) 
  // |-- [    2] --| S2(i,j,k,a) += (   -1.00000000) W218(j,a2,v0,a) W219(i,k,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W219aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no104_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO104_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W219aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W218aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no104_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO104_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W218aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no104_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO104_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W218aa_sigma_ooov_ooov.cptr(), W219aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  105] -- 
  // |-- [    0] --| W221(a3,i,a2,a) += (    1.00000000) V2(a,v0,a4,a3) W220(i,a4,a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(a3,k,j,a2) W221(a3,i,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W221aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no105_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO105_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W220aaav_sigma_ooov_ooov.cptr(), W221aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no105_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO105_X1_TYPE1_ERI_V)
    (sa, ia, W221aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  106] -- 
  // |-- [    0] --| W223(a4,a3,i,a) += (    1.00000000) V2(a,v0,a4,a3) W222(i,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(a3,a4,j,k) W223(a4,a3,i,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W223aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no106_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO106_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W222av_sigma_ooov_ooov.cptr(), W223aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no106_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO106_X1_TYPE1_ERI_V)
    (sa, ia, W223aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  107] -- 
  // |-- [    0] --| W225(a4,j,a0,a) += (    1.00000000) V2(a,v0,a4,a3) W224(a3,j,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D2(a0,a4,i,k) W225(a4,j,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W225aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no107_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO107_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W224aaav_sigma_ooov_ooov.cptr(), W225aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no107_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO107_X1_TYPE1_ERI_V)
    (sa, ia, W225aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  108] -- 
  // |-- [    0] --| W227(i,k,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) D2(a0,a2,i,k) 
  // |-- [    1] --| W226(j,a1,v0,a) += (    1.00000000) V2(v0,a,a4,a3) D2(a3,a4,j,a1) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.50000000) W226(j,a1,v0,a) W227(i,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W227aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no108_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO108_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W227aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W226aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no108_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO108_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W226aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no108_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO108_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W226aa_sigma_ooov_ooov.cptr(), W227aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  109] -- 
  // |-- [    0] --| W229(a3,i,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W228(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.33333333) D2(a3,a1,j,k) W229(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W229aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no109_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO109_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W228aaav_sigma_ooov_ooov.cptr(), W229aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no109_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO109_X1_TYPE1_ERI_V)
    (sa, ia, W229aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  110] -- 
  // |-- [    0] --| W231(a3,i,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W230(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.16666667) D2(a3,k,j,a1) W231(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W231aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no110_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO110_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W230aaav_sigma_ooov_ooov.cptr(), W231aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no110_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO110_X1_TYPE1_ERI_V)
    (sa, ia, W231aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  111] -- 
  // |-- [    0] --| W233(a3,i,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W232(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.16666667) D2(a3,a1,j,k) W233(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W233aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no111_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO111_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W232aaav_sigma_ooov_ooov.cptr(), W233aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no111_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO111_X1_TYPE1_ERI_V)
    (sa, ia, W233aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  112] -- 
  // |-- [    0] --| W235(a3,i,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W234(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.33333333) D2(a3,k,j,a1) W235(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W235aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no112_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO112_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W234aaav_sigma_ooov_ooov.cptr(), W235aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no112_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO112_X1_TYPE1_ERI_V)
    (sa, ia, W235aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  113] -- 
  // |-- [    0] --| W237(a4,a3,j,a) += (    1.00000000) V2(a,v0,a4,a3) W236(j,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    2.00000000) C2(a3,a4,i,k) W237(a4,a3,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W237aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no113_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO113_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W236av_sigma_ooov_ooov.cptr(), W237aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no113_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO113_X1_TYPE1_ERI_V)
    (sa, ia, W237aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  114] -- 
  // |-- [    0] --| W239(a3,j,a2,a) += (    1.00000000) V2(a,v0,a4,a3) W238(j,a4,a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a2,a3,k,i) W239(a3,j,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W239aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no114_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO114_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W238aaav_sigma_ooov_ooov.cptr(), W239aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no114_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO114_X1_TYPE1_ERI_V)
    (sa, ia, W239aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  115] -- 
  // |-- [    0] --| W241(j,k,a2,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,a1,j,k) 
  // |-- [    1] --| W240(i,a2,v0,a) += (    1.00000000) V2(v0,a,a4,a3) C2(a2,i,a4,a3) 
  // |-- [    2] --| S2(i,j,k,a) += (   -1.00000000) W240(i,a2,v0,a) W241(j,k,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W241aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no115_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO115_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W241aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W240aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no115_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO115_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W240aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no115_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO115_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W240aa_sigma_ooov_ooov.cptr(), W241aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  116] -- 
  // |-- [    0] --| W243(a3,j,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W242(j,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a1,a3,k,i) W243(a3,j,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W243aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no116_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO116_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W242aaav_sigma_ooov_ooov.cptr(), W243aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no116_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO116_X1_TYPE1_ERI_V)
    (sa, ia, W243aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  117] -- 
  // |-- [    0] --| W245(j,k,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,k,j,a2) 
  // |-- [    1] --| W244(i,a1,v0,a) += (    1.00000000) V2(v0,a,a4,a3) C2(a1,i,a4,a3) 
  // |-- [    2] --| S2(i,j,k,a) += (   -1.00000000) W244(i,a1,v0,a) W245(j,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W245aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no117_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO117_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W245aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W244aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no117_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO117_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W244aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no117_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO117_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W244aa_sigma_ooov_ooov.cptr(), W245aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  118] -- 
  // |-- [    0] --| W247(a4,i,a0,a) += (    1.00000000) V2(a,v0,a4,a3) W246(a3,i,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.66666667) C2(a0,a4,j,k) W247(a4,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W247aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no118_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO118_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W246aaav_sigma_ooov_ooov.cptr(), W247aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no118_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO118_X1_TYPE1_ERI_V)
    (sa, ia, W247aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  119] -- 
  // |-- [    0] --| W249(a4,i,a0,a) += (    1.00000000) V2(a,v0,a4,a3) W248(a3,i,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.33333333) C2(a0,a4,j,k) W249(a4,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W249aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no119_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO119_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W248aaav_sigma_ooov_ooov.cptr(), W249aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no119_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO119_X1_TYPE1_ERI_V)
    (sa, ia, W249aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  120] -- 
  // |-- [    0] --| W251(a4,i,a0,a) += (    1.00000000) V2(a,v0,a4,a3) W250(a3,i,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.33333333) C2(a0,k,j,a4) W251(a4,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W251aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no120_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO120_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W250aaav_sigma_ooov_ooov.cptr(), W251aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no120_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO120_X1_TYPE1_ERI_V)
    (sa, ia, W251aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  121] -- 
  // |-- [    0] --| W253(a4,i,a0,a) += (    1.00000000) V2(a,v0,a4,a3) W252(a3,i,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.66666667) C2(a0,k,j,a4) W253(a4,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W253aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no121_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO121_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W252aaav_sigma_ooov_ooov.cptr(), W253aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no121_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO121_X1_TYPE1_ERI_V)
    (sa, ia, W253aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  122] -- 
  // |-- [    0] --| W254(a3,a4,a2,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,a1,a3,a4) 
  // |-- [    1] --| W255(a2,a) += (    1.00000000) V2(v0,a,a4,a3) W254(a3,a4,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W254aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no122_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO122_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W254aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no122_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO122_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W254aaa_sigma_ooov_ooov.cptr(), W255av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  123] -- 
  // |-- [    0] --| W257(a4,a) += (    1.00000000) V2(a,v0,a4,a3) W256(a3,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a4,j,k,i) W257(a4,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W257a_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no123_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO123_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W256av_sigma_ooov_ooov.cptr(), W257a_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no123_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO123_X1_TYPE1_ERI_V)
    (sa, ia, W257a_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  124] -- 
  // |-- [    0] --| W259(a4,k,a2,a) += (    1.00000000) V2(a,v0,a4,a3) W258(a3,k,a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a2,j,a4,i) W259(a4,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W259aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no124_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO124_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W258aaav_sigma_ooov_ooov.cptr(), W259aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no124_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO124_X1_TYPE1_ERI_V)
    (sa, ia, W259aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  125] -- 
  // |-- [    0] --| W260(a3,a4,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,a2,a3,a4) 
  // |-- [    1] --| W261(a1,a) += (    1.00000000) V2(v0,a,a4,a3) W260(a3,a4,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W260aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no125_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO125_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W260aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no125_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO125_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W260aaa_sigma_ooov_ooov.cptr(), W261av_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  126] -- 
  // |-- [    0] --| W263(j,i,a0,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a1,i,a2,j) 
  // |-- [    1] --| W262(a0,k,v0,a) += (    1.00000000) V2(v0,a,a4,a3) C2(a0,k,a3,a4) 
  // |-- [    2] --| S2(i,j,k,a) += (   -1.00000000) W262(a0,k,v0,a) W263(j,i,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W263aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no126_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO126_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W263aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W262aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no126_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO126_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W262aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no126_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO126_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W262aa_sigma_ooov_ooov.cptr(), W263aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  127] -- 
  // |-- [    0] --| W265(a4,k,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W264(a3,k,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.66666667) C2(a1,j,a4,i) W265(a4,k,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W265aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no127_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO127_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W264aaav_sigma_ooov_ooov.cptr(), W265aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no127_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO127_X1_TYPE1_ERI_V)
    (sa, ia, W265aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  128] -- 
  // |-- [    0] --| W267(a4,k,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W266(a3,k,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.33333333) C2(a1,i,a4,j) W267(a4,k,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W267aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no128_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO128_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W266aaav_sigma_ooov_ooov.cptr(), W267aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no128_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO128_X1_TYPE1_ERI_V)
    (sa, ia, W267aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  129] -- 
  // |-- [    0] --| W269(a4,k,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W268(a3,k,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.33333333) C2(a1,j,a4,i) W269(a4,k,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W269aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no129_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO129_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W268aaav_sigma_ooov_ooov.cptr(), W269aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no129_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO129_X1_TYPE1_ERI_V)
    (sa, ia, W269aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  130] -- 
  // |-- [    0] --| W271(a4,k,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W270(a3,k,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.66666667) C2(a1,i,a4,j) W271(a4,k,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W271aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no130_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO130_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W270aaav_sigma_ooov_ooov.cptr(), W271aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no130_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO130_X1_TYPE1_ERI_V)
    (sa, ia, W271aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  131] -- 
  // |-- [    0] --| W273(i,k,a2,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,a1,i,k) 
  // |-- [    1] --| W272(j,a2,v0,a) += (    1.00000000) V2(v0,a,a4,a3) C2(a2,j,a4,a3) 
  // |-- [    2] --| S2(i,j,k,a) += (    2.00000000) W272(j,a2,v0,a) W273(i,k,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W273aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no131_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO131_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W273aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W272aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no131_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO131_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W272aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no131_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO131_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W272aa_sigma_ooov_ooov.cptr(), W273aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  132] -- 
  // |-- [    0] --| W275(a3,i,a2,a) += (    1.00000000) V2(a,v0,a4,a3) W274(i,a4,a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a2,j,k,a3) W275(a3,i,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W275aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no132_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO132_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W274aaav_sigma_ooov_ooov.cptr(), W275aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no132_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO132_X1_TYPE1_ERI_V)
    (sa, ia, W275aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  133] -- 
  // |-- [    0] --| W277(a4,a3,i,a) += (    1.00000000) V2(a,v0,a4,a3) W276(i,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a3,a4,j,k) W277(a4,a3,i,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W277aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no133_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO133_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W276av_sigma_ooov_ooov.cptr(), W277aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no133_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO133_X1_TYPE1_ERI_V)
    (sa, ia, W277aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  134] -- 
  // |-- [    0] --| W279(a4,j,a0,a) += (    1.00000000) V2(a,v0,a4,a3) W278(a3,j,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) C2(a0,a4,i,k) W279(a4,j,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W279aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no134_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO134_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W278aaav_sigma_ooov_ooov.cptr(), W279aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no134_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO134_X1_TYPE1_ERI_V)
    (sa, ia, W279aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  135] -- 
  // |-- [    0] --| W281(i,k,a1,v0) += (    1.00000000) T2(a1,a2,a0,v0) C2(a0,a2,i,k) 
  // |-- [    1] --| W280(j,a1,v0,a) += (    1.00000000) V2(v0,a,a4,a3) C2(a1,j,a4,a3) 
  // |-- [    2] --| S2(i,j,k,a) += (   -1.00000000) W280(j,a1,v0,a) W281(i,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W281aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no135_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO135_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W281aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W280aa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no135_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO135_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W280aa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no135_x2_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO135_X2_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, W280aa_sigma_ooov_ooov.cptr(), W281aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  136] -- 
  // |-- [    0] --| W283(a3,i,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W282(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.66666667) C2(a1,a3,k,j) W283(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W283aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no136_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO136_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W282aaav_sigma_ooov_ooov.cptr(), W283aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no136_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO136_X1_TYPE1_ERI_V)
    (sa, ia, W283aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  137] -- 
  // |-- [    0] --| W285(a3,i,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W284(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.33333333) C2(a1,j,k,a3) W285(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W285aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no137_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO137_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W284aaav_sigma_ooov_ooov.cptr(), W285aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no137_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO137_X1_TYPE1_ERI_V)
    (sa, ia, W285aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  138] -- 
  // |-- [    0] --| W287(a3,i,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W286(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.33333333) C2(a1,a3,k,j) W287(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W287aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no138_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO138_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W286aaav_sigma_ooov_ooov.cptr(), W287aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no138_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO138_X1_TYPE1_ERI_V)
    (sa, ia, W287aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  139] -- 
  // |-- [    0] --| W289(a3,i,a1,a) += (    1.00000000) V2(a,v0,a4,a3) W288(i,a4,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.66666667) C2(a1,j,k,a3) W289(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W289aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_ooov_no139_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO139_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W288aaav_sigma_ooov_ooov.cptr(), W289aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_ooov_no139_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO139_X1_TYPE1_ERI_V)
    (sa, ia, W289aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  140] -- 
  // |-- [    0] --| W290(j,i,a3,v0) += (    1.00000000) T2(a0,a1,a2,v0) D3(j,a1,i,a3,a2,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) V2(v0,a,k,a3) W290(j,i,a3,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W290aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no140_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO140_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W290aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no140_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO140_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W290aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  141] -- 
  // |-- [    0] --| W291(a1,a0,k,a3,a2,a) += (    1.00000000) V2(a,a3,v0,a2) T2(a1,a0,k,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) D3(j,a3,i,a1,a2,a0) W291(a1,a0,k,a3,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W291aaaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ooov_ooov_no141_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO141_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W291aaaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ooov_ooov_no141_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO141_X1_TYPE1_ERI_V)
    (sa, ia, W291aaaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  142] -- 
  // |-- [    0] --| W292(a1,a0,k,a3,a2,a) += (    1.00000000) V2(a,v0,a3,a2) T2(a1,a0,v0,k) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) D3(j,a1,i,a0,a3,a2) W292(a1,a0,k,a3,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sk = 0;sk < nir;++sk){ 
  for(int ik = symblockinfo.psym()(sk,I_O,I_BEGIN);ik <= symblockinfo.psym()(sk,I_O,I_END);++ik){ 
    T2b = T2.get_amp2(ik);
    orz::DTensor W292aaaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sk^sa));
    FC_FUNC(g_if_sigma_ooov_ooov_no142_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO142_X0_TYPE1_ERI_V)
      (sa, ia, sk, ik, T2b.cptr(), V2_sym.cptr(), W292aaaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_ooov_no142_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO142_X1_TYPE1_ERI_V)
      (sa, ia, sk, ik, W292aaaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[k, "active"] [notNeeded]
  } // End ik
  } // End sk
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  143] -- 
  // |-- [    0] --| W293(j,i,a3,v0) += (    1.00000000) T2(a0,a1,a2,v0) D3(j,a3,i,a1,a2,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) V2(v0,k,a3,a) W293(j,i,a3,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W293aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ooov_ooov_no143_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO143_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W293aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no143_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO143_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W293aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  144] -- 
  // |-- [    0] --| W294(a1,a0,a2,a) += (    1.00000000) V2(a,a3,v0,a2) T2(a1,a0,a3,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) D3(k,i,a1,j,a0,a2) W294(a1,a0,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W294aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ooov_ooov_no144_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO144_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W294aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ooov_ooov_no144_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO144_X1_TYPE1_ERI_V)
    (sa, ia, W294aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  145] -- 
  // |-- [    0] --| W295(a1,a0,a3,a) += (    1.00000000) V2(a,v0,a3,a2) T2(a1,a0,v0,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) D3(k,i,a1,j,a0,a3) W295(a1,a0,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W295aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_ooov_ooov_no145_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO145_X0_TYPE1_ERI_V)
      (sa, ia, sa2, ia2, T2b.cptr(), V2_sym.cptr(), W295aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_ooov_ooov_no145_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO145_X1_TYPE1_ERI_V)
    (sa, ia, W295aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  146] -- 
  // |-- [    0] --| W296(a1,a0,k,a) += (    1.00000000) V2(a,v0,k,a2) T2(a1,a0,v0,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) D2(j,a1,i,a0) W296(a1,a0,k,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W296aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_ooov_ooov_no146_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO146_X0_TYPE1_ERI_V)
      (sa, ia, sa2, ia2, T2b.cptr(), V2_sym.cptr(), W296aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_ooov_ooov_no146_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO146_X1_TYPE1_ERI_V)
    (sa, ia, W296aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  147] -- 
  // |-- [    0] --| W297(a1,a0,k,a) += (    1.00000000) V2(a,a2,v0,k) T2(a1,a0,a2,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) D2(j,a1,i,a0) W297(a1,a0,k,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W297aaa_sigma_ooov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ooov_ooov_no147_x0_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO147_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W297aaa_sigma_ooov_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ooov_ooov_no147_x1_type1_eri_v,G_IF_SIGMA_OOOV_OOOV_NO147_X1_TYPE1_ERI_V)
    (sa, ia, W297aaa_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

//-@type(1).contraction(end)

  } // End myrank
  }
  }
  orz::world().barrier();

//-@type(2).contraction(begin)
  // -- Title : sigma_ooov_ooov
  //*-- Entering to take the type 2 contractions --*//

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(i,j,k,a) += (    1.00000000) D1(i,k) W11(j,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x0_type2_eri_v,G_IF_SIGMA_OOOV_OOOV_NO0_X0_TYPE2_ERI_V)
      (sa, ia, W11av_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(i,j,k,a) += (   -0.50000000) D1(j,k) W25(i,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x0_type2_eri_v,G_IF_SIGMA_OOOV_OOOV_NO1_X0_TYPE2_ERI_V)
      (sa, ia, W25av_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(i,j,k,a) += (    0.50000000) D2(j,a2,i,k) W45(a2,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no2_x0_type2_eri_v,G_IF_SIGMA_OOOV_OOOV_NO2_X0_TYPE2_ERI_V)
      (sa, ia, W45av_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(i,j,k,a) += (    0.50000000) D2(j,a1,i,k) W49(a1,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no3_x0_type2_eri_v,G_IF_SIGMA_OOOV_OOOV_NO3_X0_TYPE2_ERI_V)
      (sa, ia, W49av_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(i,j,k,a) += (   -1.00000000) C2(a2,j,k,i) W99(a2,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no4_x0_type2_eri_v,G_IF_SIGMA_OOOV_OOOV_NO4_X0_TYPE2_ERI_V)
      (sa, ia, W99av_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(i,j,k,a) += (   -1.00000000) C2(a1,j,k,i) W103(a1,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no5_x0_type2_eri_v,G_IF_SIGMA_OOOV_OOOV_NO5_X0_TYPE2_ERI_V)
      (sa, ia, W103av_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| S2(i,j,k,a) += (    1.00000000) D1(i,k) W151(j,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no6_x0_type2_eri_v,G_IF_SIGMA_OOOV_OOOV_NO6_X0_TYPE2_ERI_V)
      (sa, ia, W151av_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| S2(i,j,k,a) += (   -0.50000000) D1(j,k) W173(i,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no7_x0_type2_eri_v,G_IF_SIGMA_OOOV_OOOV_NO7_X0_TYPE2_ERI_V)
      (sa, ia, W173av_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(i,j,k,a) += (   -1.00000000) D2(j,a2,i,k) W201(a2,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no8_x0_type2_eri_v,G_IF_SIGMA_OOOV_OOOV_NO8_X0_TYPE2_ERI_V)
      (sa, ia, W201av_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(i,j,k,a) += (    0.50000000) D2(j,a1,i,k) W207(a1,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no9_x0_type2_eri_v,G_IF_SIGMA_OOOV_OOOV_NO9_X0_TYPE2_ERI_V)
      (sa, ia, W207av_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(i,j,k,a) += (    2.00000000) C2(a2,j,k,i) W255(a2,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no10_x0_type2_eri_v,G_IF_SIGMA_OOOV_OOOV_NO10_X0_TYPE2_ERI_V)
      (sa, ia, W255av_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(i,j,k,a) += (   -1.00000000) C2(a1,j,k,i) W261(a1,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no11_x0_type2_eri_v,G_IF_SIGMA_OOOV_OOOV_NO11_X0_TYPE2_ERI_V)
      (sa, ia, W261av_sigma_ooov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
//-@type(2).contraction(end)
  orz::world().barrier();

  } // End Femto
  //*-- FEMTO ends --//*

//-@loadERI(v,end)

//-@ERI.contractions(end)

//-@D4C.contractions(begin)

//-@loadD4C(a,begin)
  //*-- FEMTO begins --//*
  // Label : d4c_o
  {


  for(int s_eri = 0;s_eri < nir;++s_eri){ 
  for(int i_eri = symblockinfo.psym()(s_eri,I_O,I_BEGIN);i_eri <= symblockinfo.psym()(s_eri,I_O,I_END);++i_eri){ 
  if(hintmo.iproc_havingimo()[i_eri] == myrank) {           
  orz::DTensor C5;
  orz::LoadBin(ctinp.dir()/(format("D4C_g[%d]")%i_eri).str()) >> C5;

  //*-- Entering to take the type 1 contractions --*//
//-@type(1).contraction(begin)
  // -- Title : sigma_ooov_ooov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(i,j,k,a) += (    0.50000000) T2(a0,a1,a2,a) C5(k,i,a1,j,a0,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x0_type1_d4c_o,G_IF_SIGMA_OOOV_OOOV_NO0_X0_TYPE1_D4C_O)
      (sa, ia, sa2, ia2, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(i,j,k,a) += (   -0.50000000) T2(a2,a1,a0,a) C5(j,a1,i,k,a0,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x0_type1_d4c_o,G_IF_SIGMA_OOOV_OOOV_NO1_X0_TYPE1_D4C_O)
      (sa, ia, sa2, ia2, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(i,j,k,a) += (   -0.50000000) T2(a1,a3,a0,a) C5(a0,a1,i,k,j,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no2_x0_type1_d4c_o,G_IF_SIGMA_OOOV_OOOV_NO2_X0_TYPE1_D4C_O)
      (sa, ia, sa3, ia3, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(i,j,k,a) += (    0.50000000) T2(a0,a1,a2,a) C5(a2,a0,j,a1,i,k) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no3_x0_type1_d4c_o,G_IF_SIGMA_OOOV_OOOV_NO3_X0_TYPE1_D4C_O)
      (sa, ia, sk, ik, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(i,j,k,a) += (   -0.50000000) T2(a1,a2,a0,a) C5(a1,a0,a2,j,k,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no4_x0_type1_d4c_o,G_IF_SIGMA_OOOV_OOOV_NO4_X0_TYPE1_D4C_O)
      (sa, ia, si, ii, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(i,j,k,a) += (   -0.50000000) T2(a1,a2,a0,a) C5(a1,a0,k,i,a2,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_ooov_no5_x0_type1_d4c_o,G_IF_SIGMA_OOOV_OOOV_NO5_X0_TYPE1_D4C_O)
      (sa, ia, sj, ij, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@loadD4C(a,end)

//-@D4C.contractions(end)

  // Do timing!
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ooov_ooov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
